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

from alchemy import AbsoluteAlchemicalFactory, AlchemicalState
import blues.ncmc as blues
import blues.ncmc_switching as blues_switching

# For parallel, import and inherit from ParallelOEMolComputeCube
class OpenMMComplexSetup(OEMolComputeCube):
    title = "Set up complex for simulation in OpenMM"
    description = """
    Set up protein:ligand complex for simulation with OpenMM.

    This cube will generate an OpenMM System object containing
    a TIP3P solvated protein:ligand complex. The ligand will be parameterized
    the smirff99Frosst.ffxml parameters, which is parsed with smarty. The complex
    will be stored into a complex.oeb.gz file and streamed into the OpenMMSimulation cube.
    """
    classification = [
        ["OpenMM", "ProtLigComplex Setup"],
    ]
    tags = [tag for lists in classification for tag in lists]

    protein = parameter.DataSetInputParameter(
        'protein',
        required=True,
        help_text='Protein PDB file')

    pH = parameter.DecimalParameter(
        'pH',
        default=7.0,
        help_text="Solvent pH used to select appropriate protein protonation state.",
    )

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)",
    )

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50,
        help_text="Salt concentration (millimolar)",
    )

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
    #    required=True,
        help_text='Forcefield parameters for molecule'
    )

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
    #    required=True,
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein'
    )

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
    #   required=True,
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent'
    )

    def begin(self):
        pdbfilename = 'protein.pdb'

        # Write the protein to a PDB
        if in_orion():
            stream = StreamingDataset(self.args.protein, input_format=".pdb")
            stream.download_to_file(pdbfilename)
        else:
            protein = oechem.OEMol()
            with oechem.oemolistream(self.args.protein) as ifs:
                if not oechem.OEReadMolecule(ifs, protein):
                    raise RuntimeError("Error reading molecule")
                with oechem.oemolostream(pdbfilename) as ofs:
                    res = oechem.OEWriteConstMolecule(ofs, protein)
                    if res != oechem.OEWriteMolReturnCode_Success:
                        raise RuntimeError("Error writing protein: {}".format(res))

        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            cubename = '[{}]'.format( str(self.name) )
            self.log.info(cubename+'Parameterizing the ligand...')
            # Generate smarty ligand structure
            from smarty.forcefield import ForceField
            ffxml = mol.GetData(oechem.OEGetTag('forcefield')).encode()
            with open('mol_parameters.ffxml', 'wb') as out:
                out.write(ffxml)
            mol_ff = ForceField(open('mol_parameters.ffxml'))
            mol_top, mol_sys, mol_pos = smarty.forcefield_utils.create_system_from_molecule(mol_ff, mol)
            molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
            #Alter molecule residue name for easy selection
            molecule_structure.residues[0].name = "MOL"
            self.log.info('\t{}'.format(str(molecule_structure)))

            self.log.info(cubename+'Parameterizing the protein...')
            #Generate protein Structure object
            forcefield = app.ForceField(self.args.protein_forcefield, self.args.solvent_forcefield)
            protein_system = forcefield.createSystem( self.proteinpdb.topology )
            protein_structure = parmed.openmm.load_topology( self.proteinpdb.topology,
                                                             protein_system,
                                                             xyz=self.proteinpdb.positions )

            self.log.info('\t{}'.format(str(protein_structure)))
            # Merge structures
            pl_structure = protein_structure + molecule_structure

            #Concatenate positions arrays (ensures same units)
            positions_unit = unit.angstroms
            positions0_dimensionless = np.array( self.proteinpdb.positions / positions_unit )
            positions1_dimensionless = np.array( molecule_structure.positions / positions_unit )
            coordinates = np.vstack((positions0_dimensionless,positions1_dimensionless))
            natoms = len(coordinates)
            positions = np.zeros([natoms,3], np.float32)
            for index in range(natoms):
                (x,y,z) = coordinates[index]
                positions[index,0] = x
                positions[index,1] = y
                positions[index,2] = z
            positions = unit.Quantity(positions, positions_unit)

            #Store Structure object
            pl_structure.positions = positions
            pl_structure.save('pl_tmp.pdb', overwrite=True)

            # Solvate with PDBFixer
            self.log.info(cubename+'Solvating system with PDBFixer...')
            self.log.info('\tpH = {}'.format(self.args.pH))
            self.log.info('\tpadding = {}'.format(unit.Quantity(self.args.solvent_padding, unit.angstroms)))
            self.log.info('\tionicStrength = {}'.format(unit.Quantity(self.args.salt_concentration, unit.millimolar)))
            fixer = pdbfixer.PDBFixer('pl_tmp.pdb')
            #fixer.findMissingResidues()
            #fixer.findMissingAtoms()
            #fixer.addMissingAtoms()
            fixer.addMissingHydrogens(self.args.pH)
            fixer.addSolvent(padding=unit.Quantity(self.args.solvent_padding, unit.angstroms),
                            ionicStrength=unit.Quantity(self.args.salt_concentration, unit.millimolar)
                            )

            # Load PDBFixer object back to Structure
            tmp = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)
            #Store positions, topology, and box vectors for solvated system
            full_positions = tmp.positions
            full_topology = tmp.topology
            full_box = tmp.box

            # Remove ligand from protein Structure by AmberMask selection
            tmp.strip(":MOL")
            tmp.save('nomol_tmp.pdb', overwrite=True)
            nomol = parmed.load_file('nomol_tmp.pdb')

            # Regenerate openMM System to parameterize solvent
            nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)

            # Regenerate parameterized solvated protein structure
            solv_structure = parmed.openmm.load_topology( nomol.topology,
                                                            nomol_system,
                                                            xyz=nomol.positions,
                                                            box=nomol.box )

            # Remerge with ligand structure
            full_structure = solv_structure + molecule_structure
            # Restore box dimensions
            full_structure.box = nomol.box
            # Save full structure
            full_structure.save('complex.pdb', overwrite=True)
            self.log.info(cubename+'Saving the protien:ligand complex')
            self.log.info('\t{}'.format(str(full_structure)))
            self.log.info('\tBox = {}'.format(full_structure.box))

            self.log.info(cubename+'Generating the OpenMM system...')
            # Regenerate OpenMM system with parmed
            system = full_structure.createSystem(nonbondedMethod=app.PME,
                                                nonbondedCutoff=10.0*unit.angstroms,
                                                constraints=app.HBonds)

            self.log.info(cubename+'Saving System to complex.oeb.gz')
            # Pack complex into oeb and emit system
            complex_mol = oechem.OEMol()
            output = OpenMMSystemOutput('output')
            with oechem.oemolistream('complex.pdb') as ifs:
                if not oechem.OEReadMolecule(ifs, complex_mol):
                    raise RuntimeError("Error reading complex pdb")
                with oechem.oemolostream('complex.oeb.gz') as ofs:
                    complex_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                    oechem.OEWriteConstMolecule(ofs, complex_mol)
            self.success.emit(complex_mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

    def end(self):
        #Clean up
        os.remove('protein.pdb')
        os.remove('mol_parameters.ffxml')
        os.remove('pl_tmp.pdb')
        os.remove('nomol_tmp.pdb')
        os.remove('complex.pdb')


class OpenMMSimulation(OEMolComputeCube):
    title = "Run simulation in OpenMM"
    description = """
    Run simulation with OpenMM for protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the protein:ligand complex, reconstruct the OpenMM System object,
    minimize the system for a max of 20 iterations (for faster runtime),
    and run 1000 MD steps at 300K. The potential energies are evaluated every
    100 steps and stored to a log file.
    The OpenMM System, State, and log file are attached to the OEMol and saved
    to the file simulation.oeb.gz.

    The simulation.oeb.gz file, containing the State can then be reused to
    restart the MD simulation.
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

    steps = parameter.IntegerParameter(
        'steps',
        default=1000,
        help_text="Number of MD steps")

    complex_mol = parameter.DataSetInputParameter(
        'complex_mol',
        default='complex.oeb.gz',
        #required=True,
        help_text='Single protein to Dock Against')

    def begin(self):
        # Initialize openmm integrator
        self.integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)

    def process(self, mol, port):
        try:
            cubename = '[{}]'.format( str(self.name) )
            # Reconstruct the OpenMM system
            if 'system' in mol.GetData().keys():
                self.log.info(cubename+'Regenerating System from OEMol')
                serialized_system = mol.GetData(oechem.OEGetTag('system'))
                system = openmm.XmlSerializer.deserialize( serialized_system )

                self.log.info(cubename+'Regenerating positions and topology from OEMol')
                positions = unit.Quantity(np.zeros([mol.NumAtoms(), 3], np.float32), unit.angstroms)
                coords = mol.GetCoords()
                for index in range(mol.NumAtoms()):
                    positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
                topology = forcefield_generators.generateTopologyFromOEMol(mol)
                self.log.info(str(topology))

                # Initialize Simulation
                simulation = app.Simulation(topology, system, self.integrator)
                #simulation = app.Simulation(topology, system, self.integrator, openmm.Platform.getPlatformByName('CPU'))
                platform = simulation.context.getPlatform().getName()
                self.log.info(cubename+'Running OpenMMSimulation on Platform {}'.format(platform))
            else:
                raise RuntimeError('Could not find system from mol')

            # Check if mol has State data attached
            if 'state' in mol.GetData().keys():
                self.log.info(cubename+'Found a saved State, restarting simulation')
                outfname = 'restart'
                mol.GetData(oechem.OEGetTag('state'))
                serialized_state = mol.GetData(oechem.OEGetTag('state'))
                state = openmm.XmlSerializer.deserialize( serialized_state )
                simulation.context.setState(state)
            else:
                self.log.info(cubename+'Minimizing system...')
                outfname = 'simulation'
                # Set initial positions and velocities then minimize
                simulation.context.setPositions(positions)
                simulation.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)

                # Temporarily, place some restrictions on minization to run faster
                if in_orion():
                    simulation.minimizeEnergy(tolerance=unit.Quantity(10.0,unit.kilojoules/unit.moles),maxIterations=100)
                else:
                    simulation.minimizeEnergy()
                st = simulation.context.getState(getPositions=True,getEnergy=True)
                self.log.info('\tMinimized energy is {}'.format(st.getPotentialEnergy()))

            self.log.info(cubename+'Running {} MD steps at {}K'.format(self.args.steps, self.args.temperature))
            # Do MD simulation and report energies
            statereporter = app.StateDataReporter(outfname+'.log', 100, step=True, potentialEnergy=True, temperature=True)
            simulation.reporters.append(statereporter)
            outlog = open(outfname+'.log', 'r')
            simulation.step(self.args.steps)
            self.log.info(outlog.read())

            # Save serialized State object
            state = simulation.context.getState( getPositions=True,
                                              getVelocities=True,
                                              getParameters=True )

            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            self.log.info(cubename+'Saving to {}'.format(outfname+'.oeb.gz'))
            with oechem.oemolostream(outfname+'.oeb.gz') as ofs:
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
