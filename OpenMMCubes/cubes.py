import io, os, time, traceback, base64, smarty, parmed, pdbfixer, mdtraj, pickle
from openeye import oechem
from sys import stdout
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app

from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort, ParallelOEMolComputeCube
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
import OpenMMCubes.utils as utils
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename


class OpenMMComplexSetup(OEMolComputeCube):
    title = "OpenMMComplexSetup"
    description = """
    Set up protein:ligand complex for simulation with OpenMM.

    This cube will generate an OpenMM System containing
    a TIP3P solvated protein:ligand complex. The complex
    will be stored into a <idtag>-complex.oeb.gz file, with the System and Structure
    attached and streamed into the OpenMMSimulation cube.
    """
    classification = ["Complex Setup"]
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

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

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein'
    )

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent'
    )

    def begin(self):
        pdbfilename = 'protein.pdb'
        protein = oechem.OEMol()
        self.args.protein = download_dataset_to_file(self.args.protein)
        with oechem.oemolistream(self.args.protein) as ifs:
            if not oechem.OEReadMolecule(ifs, protein):
                raise RuntimeError("Error reading protein")
        with oechem.oemolostream(pdbfilename) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, protein)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing protein: {}".format(res))

        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            if utils.PackageOEMol.checkTags(mol, ['idtag', 'structure']):
                data = utils.PackageOEMol.unpack(mol)
                outfname = '{}-complex'.format(data['idtag'])

            #Generate protein Structure object
            forcefield = app.ForceField(self.args.protein_forcefield, self.args.solvent_forcefield)
            protein_system = forcefield.createSystem( self.proteinpdb.topology )
            protein_structure = parmed.openmm.load_topology( self.proteinpdb.topology,
                                                             protein_system,
                                                             xyz=self.proteinpdb.positions )

            # Merge structures to prevent adding solvent in pocket
            pl_structure = protein_structure + data['structure']
            self.log.info('{}-complex: {}'.format(data['idtag'], pl_structure))

            # Retain positions and save
            pl_structure.positions = utils.combinePostions(protein_structure.positions,
                                            data['structure'].positions)
            pl_structure.save(outfname+'-pl.tmp',format='pdb',overwrite=True)

            # Solvate with PDBFixer
            self.log.info('PDBFixer solvating {}-complex:'.format(data['idtag']))
            self.log.info('\tpH = {}'.format(self.args.pH))
            self.log.info('\tpadding = {}'.format(unit.Quantity(self.args.solvent_padding, unit.angstroms)))
            self.log.info('\tionicStrength = {}'.format(unit.Quantity(self.args.salt_concentration, unit.millimolar)))
            fixer = pdbfixer.PDBFixer(outfname+'-pl.tmp')
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.replaceNonstandardResidues()
            #fixer.removeHeterogens(False)
            fixer.addMissingAtoms()
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
            tmp.strip(":LIG")
            tmp.save(outfname+'-nomol.tmp',format='pdb',overwrite=True)
            # Reload PDBFile
            nomol = app.PDBFile(outfname+'-nomol.tmp')
            nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)
            # Regenerate parameterized solvated protein structure
            solv_structure = parmed.openmm.load_topology(nomol.topology,
                                                        nomol_system,
                                                        xyz=nomol.positions,
                                                        box=full_box)

            # Remerge with ligand structure
            full_structure = solv_structure + data['structure']
            # Restore box dimensions
            full_structure.box = full_box
            # Save full structure
            full_structure.save(outfname+'.pdb', overwrite=True)
            self.log.info('Solvated {}-complex {}'.format(data['idtag'], full_structure))
            self.log.info('\tBox = {}'.format(full_structure.box))

            # Regenerate OpenMM system with parmed
            system = full_structure.createSystem(nonbondedMethod=app.PME,
                                                nonbondedCutoff=10.0*unit.angstroms,
                                                constraints=app.HBonds)

            # Pack solvated complex into oeb and emit system
            complex_mol = oechem.OEMol()
            sys_out = OpenMMSystemOutput('sys_out')
            struct_out = ParmEdStructureOutput('struct_out')

            pdbout = open(outfname+'.pdb', 'rb')
            with oechem.oemolistream(outfname+'.pdb') as ifs:
                if not oechem.OEReadMolecule(ifs, complex_mol):
                    raise RuntimeError("Error reading {}.pdb".format(outfname))
                complex_mol.SetData(oechem.OEGetTag('idtag'), data['idtag'])
                complex_mol.SetData(oechem.OEGetTag('system'), sys_out.encode(system))
                complex_mol.SetData(oechem.OEGetTag('structure'), struct_out.encode(full_structure))
                complex_mol.SetData(oechem.OEGetTag('pdb'), pdbout.read())
            self.success.emit(complex_mol)
            pdbout.close()

            #Cleanup
            tmpfiles = ['protein.pdb', outfname+'-pl.tmp',
                       outfname+'-nomol.tmp', outfname+'.pdb']
            for f in tmpfiles:
                os.remove(f)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class OpenMMSimulation(OEMolComputeCube):
    title = "Run simulation in OpenMM"
    description = """
    Run simulation with OpenMM for protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the protein:ligand complex, reconstruct the OpenMM System,
    minimize the system, save the minimized PDB, and run 1000 MD steps at 300K.
    The potential energies are evaluated every 1000 steps and stored to a log file.
    Stdout is a progress/benchmark timings reporter every 1000 steps.
    The Structure, OpenMM System, State, and log file are attached to the OEMol and
    saved to the file simulation.oeb.gz.

    The simulation.oeb.gz file, containing the State can then be reused to
    restart the MD simulation.
    """
    classification = ["Simulation"]
    tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)"
    )
    steps = parameter.IntegerParameter(
        'steps',
        default=50000,
        help_text="Number of MD steps")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=1000,
        help_text="Step interval for reporting data."
    )

    def check_tagdata(self, mol):
        if 'idtag' not in mol.GetData().keys():
            raise RuntimeError('Could not find idtag for molecule')
        else:
            idtag =  mol.GetData(oechem.OEGetTag('idtag'))
            if not os.path.exists('./output'):
                os.makedirs('./output')
            self.outfname = 'output/{}-simulation'.format(idtag)
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
            self.outfname = 'output/{}-restart'.format(self.idtag)
        else:
            self.state = None

        if not any([self.idtag, self.system, self.structure]):
            raise RuntimeError('Missing tagged generic data')
        else:
            return True

    def setReporters(self):
        progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                            reportInterval=self.args.reporter_interval,
                                            totalSteps=self.args.steps,
                                            time=True, speed=True, progress=True,
                                            elapsedTime=True, remainingTime=True)
        #progress_reporter = parmed.openmm.reporters.ProgressReporter(stdout, totalSteps=self.args.steps, reportInterval=self.args.reporter_interval)

        #state_reporter = app.StateDataReporter(self.outfname+'.log', separator="\t",
        #                                    reportInterval=self.args.reporter_interval,
        #                                    step=True,
        #                                    potentialEnergy=True, totalEnergy=True,
        #                                    volume=True, temperature=True)
        state_reporter = parmed.openmm.reporters.StateDataReporter(self.outfname+'.log', reportInterval=self.args.reporter_interval)
        traj_reporter = mdtraj.reporters.HDF5Reporter(self.outfname+'.h5', self.args.reporter_interval)
        #cdf_reporter = mdtraj.reporters.NetCDFReporter(self.outfname+'.nc', self.args.reporter_interval)
        self.reporters = [progress_reporter, state_reporter, traj_reporter]#, cdf_reporter]
        return self.reporters

    def begin(self):
        pass

    def process(self, complex_mol, port):
        try:
            outmol = oechem.OEMol()
            if self.check_tagdata(complex_mol):
                idtag = self.idtag
                outfname = self.outfname
                system = self.system
                structure = self.structure
                positions = structure.positions
                topology = structure.topology
            # Initialize Simulation
            integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
            simulation = app.Simulation(topology, system, integrator)
            platform = simulation.context.getPlatform().getName()
            self.log.info('Running OpenMMSimulation on Platform {}'.format(platform))

            # Check if mol has State data attached
            if self.state:
                simulation.context.setState(self.state)
            else:
                # Set initial positions and velocities then minimize
                simulation.context.setPositions(positions)
                simulation.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
                #init = simulation.context.getState(getEnergy=True)
                #self.log.info('Initial energy is {}'.format(init.getPotentialEnergy()))
                self.log.info('Minimizing {} system...'.format(idtag))
                min_reporter = parmed.openmm.reporters.EnergyMinimizerReporter(stdout)
                simulation.reporters.append(min_reporter)
                simulation.minimizeEnergy()
                min_state = simulation.context.getState(getPositions=True,getEnergy=True)
                self.log.info('Minimized energy is {}'.format(min_state.getPotentialEnergy()))
                minfname = 'output/{}-minimized.pdb'.format(idtag)
                with open(minfname, 'w') as minout:
                    app.PDBFile.writeFile(simulation.topology, min_state.getPositions(), minout)

                with open(minfname, 'rb') as minin:
                    outmol.SetData(oechem.OEGetTag('pdb') , minin.read())
                os.remove(minfname)

            #Append Reporters to simulation
            reporters = self.setReporters()
            for rep in reporters:
                simulation.reporters.append(rep)

            self.log.info('Running {} MD steps at {}K'.format(self.args.steps, self.args.temperature))
            simulation.step(self.args.steps)
            outlog = open(outfname+'.log', 'r')
            self.log.info(outlog.read())

            # Save serialized State object
            state = simulation.context.getState(getPositions=True,
                                                getVelocities=True,
                                                getParameters=True,
                                                getForces=True,
                                                getParameterDerivatives=True,
                                                getEnergy=True,
                                                enforcePeriodicBox=True)

            pdbfname = 'output/{}-simulation.pdb'.format(idtag)
            with open(pdbfname, 'w') as simout:
                app.PDBFile.writeFile(simulation.topology, state.getPositions(), simout)
            for rep in reporters:
                try:
                    print(rep)
                    rep.close()
                except Exception as e:
                    print(e)
                    pass
            h5traj = mdtraj.load(outfname+'.h5')
            pkl_traj = pickle.dumps(h5traj)
            b64_traj = base64.b64encode(pkl_traj)

            # Attach openmm objects to newmol, emit to output
            pdbout = open(pdbfname, 'rb')
            output = OpenMMSystemOutput('output')
            with oechem.oemolistream(pdbfname) as ifs:
                if not oechem.OEReadMolecule(ifs, outmol):
                    raise RuntimeError("Error reading PDB: {}".format(pdbfname))
                outmol.SetData(oechem.OEGetTag('idtag'), idtag)
                outmol.SetData(oechem.OEGetTag('system'), output.encode(system))
                outmol.SetData(oechem.OEGetTag('state'), output.encode(state))
                outmol.SetData(oechem.OEGetTag('pdb'), pdbout.read())
                outmol.SetData(oechem.OEGetTag('traj'), b64_traj)
                outmol.SetData(oechem.OEGetTag('log'), outlog.read())

            self.success.emit(outmol)

            outlog.close()
            pdbout.close()
            os.remove(pdbfname)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                complex_mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(complex_mol)
