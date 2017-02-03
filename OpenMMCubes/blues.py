import time
import traceback
import numpy as np
from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from OpenMMCubes.ports import OpenMMSystemOutput, OpenMMSystemInput
from simtk import unit, openmm
from simtk.openmm import app

from openeye import oechem
import os, smarty, parmed, pdbfixer
from openmoltools import forcefield_generators
from simtk.openmm import XmlSerializer

from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort

from alchemy import AbsoluteAlchemicalFactory, AlchemicalState
import blues.ncmc as blues
import blues.ncmc_switching as blues_switching

class BluesNCMC(OEMolComputeCube):
    title = "Run enhanced sampling method BLUES"
    description = """
    Some long description.
    """
    classification = [
        ["Testing", "OpenMM"],
        ["Testing", "Simulation"],
    ]
    tags = [tag for lists in classification for tag in lists]

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)"
    )

    friction = parameter.DecimalParameter(
        'friction',
        default=1,
        help_text="Friction coefficient"
    )

    timestep = parameter.DecimalParameter(
        'timestep',
        default=0.002,
        help_text='timestep'
    )

    mdsteps = parameter.IntegerParameter(
        'mdsteps',
        default=100,
        help_text="Number of MD steps")

    ncsteps = parameter.IntegerParameter(
        'ncsteps',
        default=10,
        help_text="Number of NCMC steps")

    nciter = parameter.IntegerParameter(
        'nciter',
        default=10,
        help_text="Number of NCMC iterations")

    complex_mol = parameter.DataSetInputParameter(
        'complex_mol',
        default='complex.oeb.gz',
        #required=True,
        help_text='Single protein to Dock Against')

    def begin(self):
        #Initialize integrators
        self.md_integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin,
                                                  self.args.friction/unit.picosecond,
                                                  self.args.timestep*unit.picoseconds)
        self.alch_integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin,
                                                  self.args.friction/unit.picosecond,
                                                  self.args.timestep*unit.picoseconds)

        #Defines ncmc move eqns for lambda peturbation of sterics/electrostatics
        self.functions = { 'lambda_sterics' : 'step(0.199999-lambda) + step(lambda-0.2)*step(0.8-lambda)*abs(lambda-0.5)*1/0.3 + step(lambda-0.800001)',
                    'lambda_electrostatics' : 'step(0.2-lambda)- 1/0.2*lambda*step(0.2-lambda) + 1/0.2*(lambda-0.8)*step(lambda-0.8)' }


    def process(self, mol, port):
        try:
            cubename = '[{}]'.format( str(self.name) )
            # Reconstruct the OpenMM system
            if 'system' in mol.GetData().keys():
                self.log.info(cubename+'Regenerating System from OEMol')
                serialized_system = mol.GetData(oechem.OEGetTag('system'))
                system = openmm.XmlSerializer.deserialize( serialized_system )

                self.log.info(cubename+'Regenerating positions and topology from OEMol.')
                positions = unit.Quantity(np.zeros([mol.NumAtoms(), 3], np.float32), unit.angstroms)
                coords = mol.GetCoords()
                for index in range(mol.NumAtoms()):
                    positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
                topology = forcefield_generators.generateTopologyFromOEMol(mol)
                self.log.info(str(topology))
            else:
                raise RuntimeError('Could not find system from mol')

            # Get indices of ligand atoms
            lig_atoms = []
            if not oechem.OEHasResidues(mol):
                oechem.OEPerceiveResidues(mol, oechem.OEPreserveResInfo_All)
            for atom in mol.GetAtoms():
                thisRes = oechem.OEAtomGetResidue(atom)
                resname = thisRes.GetName()
                if 'MOL' in resname:
                    lig_atoms.append(atom.GetIdx())
            self.log.info(cubename+'Selected ligand atoms: {}'.format(str(lig_atoms)))

            #Initialize MD simualtion
            md_sim = app.Simulation(topology, system, self.md_integrator)
            md_sim.context.setPositions(positions)
            md_sim.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
            platform = md_sim.context.getPlatform().getName()
            self.log.info(cubename+'Running BluesNCMC on Platform {}'.format(platform))

            #Add reporters
            statereporter = app.StateDataReporter('blues_mdsim.log', 100, step=True, potentialEnergy=True, temperature=True)
            md_sim.reporters.append(statereporter)
            md_sim.reporters.append(app.dcdreporter.DCDReporter('traj.dcd', self.args.mdsteps))

            #Initialize Alchemical Simulation
            # performs alchemical corrections
            # Reporter for NCMC moves
            alch_sim = app.Simulation(topology, system, self.alch_integrator)
            # Generate Alchemical System
            factory = AbsoluteAlchemicalFactory(system,
                                                ligand_atoms=lig_atoms,
                                                annihilate_sterics=True,
                                                annihilate_electrostatics=True)
            alch_system = factory.createPerturbedSystem()

            # Generate NC Integrator/Context
            nc_integrator = blues_switching.NCMCVVAlchemicalIntegrator(self.args.temperature*unit.kelvin,
                                                      alch_system,
                                                      self.functions,
                                                      nsteps=self.args.ncsteps,
                                                      direction='insert',
                                                      timestep=0.001*unit.picoseconds,
                                                      steps_per_propagation=1)
            nc_context = openmm.Context(alch_system, nc_integrator)

            #Initialize BLUES engine
            blues = blues.SimNCMC(temperature=self.args.temperature*unit.kelvin, residueList=lig_atoms)
            #Define NC Move
            # Rotation around the COM at some step
            # Again to maintain symmetry of ncmc move
            rot_step = (self.args.ncsteps/2) - 1
            nc_move = [[blues.rotationalMove, [rot_step]]]

            # actually run
            outlog = open('blues_mdsim.log', 'r')
            blues.get_particle_masses(system, residueList=lig_atoms)
            blues.runSim(md_sim, nc_context, nc_integrator,
                            alch_sim, movekey=nc_move,
                            niter=self.args.nciter, nstepsNC=self.args.ncsteps, nstepsMD=self.args.mdsteps,
                            alchemical_correction=True)

            self.log.info(outlog.read())
            state = blues.md_simulation.context.getState(getPositions=True,getEnergy=True)
            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            self.log.info(cubename+'Saving to {}'.format('blues.oeb.gz'))
            with oechem.oemolostream('blues.oeb.gz') as ofs:
                mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                mol.SetData(oechem.OEGetTag('state'), output.encode(state))
                mol.SetData(oechem.OEGetTag('log'), outlog.read())
                oechem.OEWriteConstMolecule(ofs, mol)
            self.success.emit(mol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(mol)
