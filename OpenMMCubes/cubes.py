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
import OpenMMCubes.simtools as simtools

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
        default=7.4,
        help_text="Solvent pH used to select appropriate protein protonation state.")

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)")

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=100,
        help_text="Salt concentration (millimolar)")

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein')

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent')

    def begin(self):
        pdbfilename = 'protein.pdb'
        protein = oechem.OEMol()
        self.args.protein = utils.download_dataset_to_file(self.args.protein)
        with oechem.oemolistream(self.args.protein) as ifs:
            if not oechem.OEReadMolecule(ifs, protein):
                raise RuntimeError("Error reading protein")
        with oechem.oemolostream(pdbfilename) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, protein)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing protein: {}".format(res))

        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(pdbfilename)
        if self.proteinpdb:
            utils.cleanup(['protein.pdb'])
        self.opt = vars(self.args)
        self.opt['logger'] = self.log

    def process(self, mol, port):
        try:
            # Check for generic data.
            req_tags = ['idtag', 'structure']
            if utils.PackageOEMol.checkTags(mol, req_tags):
                gd = utils.PackageOEMol.unpack(mol)
                self.opt['outfname'] = '{}-complex'.format(gd['idtag'])

            # Generate parameterized protein Structure
            protein_structure = simtools.genProteinStructure(self.proteinpdb,**self.opt)

            # Merge structures to prevent adding solvent in pocket
            # Ligand must be in docked position
            pl_structure = simtools.mergeStructure(protein_structure, gd['structure'] )
            self.log.info('{}: {}'.format(self.opt['outfname'], pl_structure))

            # Returns solvated system w/o ligand.
            solv_structure = simtools.solvateComplexStructure(pl_structure, **self.opt)

            # Remerge with ligand structure
            full_structure = simtools.mergeStructure(solv_structure, gd['structure'])
            self.log.info('Solvated {}: {}'.format(self.opt['outfname'], full_structure))
            self.log.info('\tBox = {}'.format(full_structure.box))

            # Emit OEMol with attached Structure
            packedmol = utils.PackageOEMol.pack(mol, full_structure)
            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class OpenMMSimulation(OEMolComputeCube):
    title = "OpenMMSimulation"
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
        help_text="Temperature (Kelvin)")

    steps = parameter.IntegerParameter(
        'steps',
        default=50000,
        help_text="Number of MD steps")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=1000,
        help_text="Step interval for reporting data.")

    trajectory_filetype = parameter.StringParameter(
        'trajectory_filetype',
        default='HDF5',
        help_text="NetCDF, DCD, HDF5. Filetype to write trajectory files")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The nonbonded cutoff in angstroms.
        This is ignored if nonbondedMethod is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['logger'] = self.log
    def process(self, mol, port):
        try:
            # Check for generic data.
            req_tags = ['idtag', 'structure']
            if utils.PackageOEMol.checkTags(mol, req_tags):
                gd = utils.PackageOEMol.unpack(mol)
                self.opt['outfname'] = '{}-simulation'.format(gd['idtag'])

            # Generate Simulation from Structure
            simulation = simtools.genSimFromStruct(gd['structure'], **self.opt)
            platform = simulation.context.getPlatform().getName()
            self.log.info('Running OpenMMSimulation on Platform {}'.format(platform))

            # Check if mol has State data attached
            if 'state' in gd.keys():
                self.log.info('Restarting from saved state...')
                simulation.context.setState(gd['state'])
            else:
                self.log.info('Minimizing system...')
                simulation = simtools.minimizeSimulation(simulation, **self.opt)

            for rep in simtools.getReporters(**self.opt):
                simulation.reporters.append(rep)

            self.log.info('Running {steps} MD steps at {temperature}K'.format(**self.opt))
            simulation.step(self.args.steps)

            packedmol = utils.PackageOEMol.pack(mol, simulation, **self.opt)
            self.success.emit(packedmol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(mol)
