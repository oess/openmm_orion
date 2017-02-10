import io, os, time, traceback, base64, smarty, parmed, pdbfixer
from openeye import oechem
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app

from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
import OpenMMCubes.utils as utils
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )

from alchemy import AbsoluteAlchemicalFactory, AlchemicalState
import blues.ncmc as ncmc
import blues.ncmc_switching as blues_switching

class BluesNCMC(OEMolComputeCube):
    title = "Run BLUES: Binding modes of Ligands Using Enhanced Sampling"
    description = """
    Run BLUES: Binding modes of Ligands Using Enhanced Sampling
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

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=1000,
        help_text="Step interval for reporting data."
    )

    mdsteps = parameter.IntegerParameter(
        'mdsteps',
        default=25000,
        help_text="Number of MD steps")

    ncsteps = parameter.IntegerParameter(
        'ncsteps',
        default=25,
        help_text="Number of NCMC steps")

    nciter = parameter.IntegerParameter(
        'nciter',
        default=10,
        help_text="Number of NCMC iterations")

    def check_tagdata(self, mol):
        if 'idtag' not in mol.GetData().keys():
            raise RuntimeError('Could not find idtag for molecule')
        else:
            idtag =  mol.GetData(oechem.OEGetTag('idtag'))
            self.idtag = idtag
        if 'system' not in mol.GetData().keys():
            raise RuntimeError("Could not find system for molecule")
        else:
            sys_in = OpenMMSystemInput('sys_in')
            sys_tag = oechem.OEGetTag('system')
            system = sys_in.decode(mol.GetData(sys_tag))
            self.system = system
        if 'structure' not in mol.GetData().keys():
            raise RuntimeError('Could not find structure for molecule')
        else:
            struct_in = ParmEdStructureInput('struct_in')
            struct_tag = oechem.OEGetTag('structure')
            structure = struct_in.decode(mol.GetData(struct_tag))
            self.structure = structure

        # Check if mol has State data attached
        if 'state' in mol.GetData().keys():
            self.log.info('Found a saved State, restarting simulation')
            mol.GetData(oechem.OEGetTag('state'))
            serialized_state = mol.GetData(oechem.OEGetTag('state'))
            state = openmm.XmlSerializer.deserialize( serialized_state )
            self.state = state
            self.outfname = 'output/{}-blues'.format(self.idtag)
        else:
            self.outfname = 'output/{}-blues_0'.format(self.idtag)
            oechem.OEThrow.Warning('Could not find a saved state, recommend equilibrating system before BLUES.')
            self.state = None

        if not any([self.idtag, self.system, self.structure]):
            raise RuntimeError('Missing tagged generic data')
        else:
            return True

    def atomIndexfromTop(self, resname, topology):
        ligand_atoms = []
        for atom in topology.atoms():
            if str(resname) in atom.residue.name:
                ligand_atoms.append(atom.index)
        self.ligand_atoms = ligand_atoms
        return self.ligand_atoms

    def setReporters(self):
        from sys import stdout
        #Add initial number of MD steps
        totalSteps = 50000 + int(self.args.mdsteps)
        progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                            reportInterval=self.args.reporter_interval,
                                            totalSteps=totalSteps,
                                            time=True, speed=True, progress=True,
                                            elapsedTime=True, remainingTime=True)

        state_reporter = app.StateDataReporter(self.outfname+'.log', separator="\t",
                                            reportInterval=self.args.reporter_interval,
                                            step=True,
                                            potentialEnergy=True, totalEnergy=True,
                                            volume=True, temperature=True)
        chk_reporter = app.checkpointreporter.CheckpointReporter(self.outfname+'.chk',
                                                                reportInterval=int(self.args.reporter_interval * 10))
        import mdtraj
        traj_reporter = mdtraj.reporters.HDF5Reporter(self.outfname+'.h5', self.args.reporter_interval)
        dcd_reporter = app.dcdreporter.DCDReporter(self.outfname+'.dcd', self.args.reporter_interval)
        self.reporters = [progress_reporter, state_reporter, traj_reporter, dcd_reporter, chk_reporter]
        return self.reporters

    def process(self, mol, port):
        try:
            if self.check_tagdata(mol):
                idtag = self.idtag
                outfname = self.outfname
                system = self.system
                structure = self.structure
                positions = structure.positions
                topology = structure.topology
                ligand_atoms = self.atomIndexfromTop('MOL', topology)
                if not ligand_atoms:
                    raise RuntimeError('Could not found resname MOL in molecule {}'.format(idtag))
                else:
                    self.log.info('Selected ligand atoms: {}'.format(str(ligand_atoms)))
            #Initialize integrators
            md_integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin,
                                                        self.args.friction/unit.picosecond,
                                                        self.args.timestep*unit.picoseconds)
            alch_integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin,
                                                        self.args.friction/unit.picosecond,
                                                        self.args.timestep*unit.picoseconds)

            #Defines ncmc move eqns for lambda peturbation of sterics/electrostatics
            functions = { 'lambda_sterics' : 'step(0.199999-lambda) + step(lambda-0.2)*step(0.8-lambda)*abs(lambda-0.5)*1/0.3 + step(lambda-0.800001)',
            'lambda_electrostatics' : 'step(0.2-lambda)- 1/0.2*lambda*step(0.2-lambda) + 1/0.2*(lambda-0.8)*step(lambda-0.8)' }

            #Initialize MD simualtion
            md_sim = app.Simulation(topology, system, md_integrator)
            if self.state:
                md_sim.context.setState(self.state)
            else:
                md_sim.context.setPositions(positions)
                md_sim.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
            platform = md_sim.context.getPlatform().getName()
            self.log.info('Running BluesNCMC on Platform {}'.format(platform))

            #Append Reporters to simulation
            reporters = self.setReporters()
            for rep in reporters:
                md_sim.reporters.append(rep)

            #Initialize Alchemical Simulation
            # performs alchemical corrections
            # Reporter for NCMC moves
            alch_sim = app.Simulation(topology, system, alch_integrator)
            # Generate Alchemical System
            factory = AbsoluteAlchemicalFactory(system, ligand_atoms,
                                                annihilate_sterics=True,
                                                annihilate_electrostatics=True)
            alch_system = factory.createPerturbedSystem()

            # Generate NC Integrator/Context
            nc_integrator = blues_switching.NCMCVVAlchemicalIntegrator(self.args.temperature*unit.kelvin,
                                                      alch_system, functions,
                                                      nsteps=self.args.ncsteps,
                                                      direction='insert',
                                                      timestep=0.001*unit.picoseconds,
                                                      steps_per_propagation=1)
            nc_context = openmm.Context(alch_system, nc_integrator)

            #Initialize BLUES engine
            blues = ncmc.SimNCMC(temperature=self.args.temperature*unit.kelvin, residueList=ligand_atoms)
            #Define NC Move
            # Rotation around the COM at some step
            # Again to maintain symmetry of ncmc move
            rot_step = (self.args.ncsteps/2) - 1
            nc_move = [[blues.rotationalMove, [rot_step]]]

            # actually run
            outlog = open(outfname+'.log', 'r')
            self.log.info('Running {} MD steps at {}K'.format(self.args.mdsteps, self.args.temperature))
            blues.get_particle_masses(system, residueList=ligand_atoms)
            blues.runSim(md_sim, nc_context, nc_integrator, alch_sim, movekey=nc_move,
                        niter=self.args.nciter, nstepsNC=self.args.ncsteps, nstepsMD=self.args.mdsteps,
                        alchemical_correction=True)
            self.log.info(outlog.read())

            state = blues.md_simulation.context.getState(getPositions=True,getEnergy=True)
            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            struct_out = ParmEdStructureOutput('struct_out')
            mol.SetData(oechem.OEGetTag('system'), output.encode(system))
            mol.SetData(oechem.OEGetTag('structure'), struct_out.encode(structure))
            mol.AddData(oechem.OEGetTag('state'), output.encode(state))
            mol.AddData(oechem.OEGetTag('log'), outlog.read())
            self.success.emit(mol)
            outlog.close()
        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(mol)
