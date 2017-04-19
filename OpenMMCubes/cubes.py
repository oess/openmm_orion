import io, os, traceback
import uuid
import numpy as np
import mdtraj, parmed
from openeye import oechem
from simtk import openmm, unit
from simtk.openmm import app
import OpenMMCubes.simtools as simtools
import OpenMMCubes.utils as utils
from floe.api import ParallelOEMolComputeCube, parameter

class OpenMMComplexSetup(ParallelOEMolComputeCube):
    title = "OpenMM Complex Setup"
    version = "0.0.1"
    classification = [["Protein Preparation", "OpenMM", "Forcefield Assignment"],
    ["Protein Preparation", "PDBFixer", "Solvate"],
    ["Protein Preparation", "PDBFixer", "Add Missing Atoms"],
    ["Protien Preparation", "PDBFixer", "Assign Protonation States"],
    ["Protein-Ligand Preparation", "ParmEd", "Generate Complex"]]
    tags = ['PDBFixer', 'OpenMM', 'ParmEd', 'Parallel Cube']
    description = """
    Using PDBFixer, add missing atoms, assign protonation state with a given pH,
    solvate the system with TIP3P, and assign forcefield parameters (default: amber99sbildn).
    Generate a parameterized parmed Structure of the solvated protein:ligand complex.

    Input:
    -------
    protein - Requires a PDB file of the protein.
    oechem.OEMol - Streamed-in charged and docked molecule with explicit hydrogens.

    Output:
    -------
    oechem.OEMol - Emits molecule with attachments:
        - SDData Tags: { Structure: str <parmed.Structure> }
        - Generic Tags: { Structure : parmed.Structure (base64-encoded) }
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

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
        #pdbfilename = 'protein.pdb'
        protein = oechem.OEMol()
        self.args.protein = utils.download_dataset_to_file(self.args.protein)
        # Read the PDB file into an OEMol
        with oechem.oemolistream(self.args.protein) as ifs:
            if not oechem.OEReadMolecule(ifs, protein):
                raise RuntimeError("Error reading protein")
        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(self.args.protein)

        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, mol, port):
        try:
            # Check for generic data.
            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                molecule_structure = gd['Structure']
                self.opt['outfname'] = '{}-complex'.format(gd['IDTag'])

            # Generate parameterized protein Structure
            protein_structure = simtools.genProteinStructure(
                self.proteinpdb, **self.opt)

            # Merge structures to prevent adding solvent in pocket
            # Ligand must be in docked position
            pl_structure = simtools.mergeStructure(
                protein_structure, molecule_structure)
            self.log.info('{}: {}'.format(self.opt['outfname'], pl_structure))

            # Returns solvated system w/o ligand.
            solv_structure = simtools.solvateComplexStructure(
                pl_structure, **self.opt)

            # Remerge with ligand structure
            full_structure = simtools.mergeStructure(
                solv_structure, molecule_structure)
            self.log.info('Solvated {}: {}'.format(
                self.opt['outfname'], full_structure))
            self.log.info('\tBox = {}'.format(full_structure.box))

            # Emit OEMol with attached Structure
            oechem.OESetSDData(mol, 'Structure', str(full_structure))
            packedmol = utils.PackageOEMol.pack(mol, full_structure)
            packedmol.SetData(oechem.OEGetTag(
                'outfname'), self.opt['outfname'])

            # Save the reference positions in OEMol
            ref_positions = full_structure.positions
            packedpos = utils.PackageOEMol.encodePyObj(ref_positions)
            packedmol.SetData(oechem.OEGetTag('OEMDDataRefPositions'), packedpos)
            
            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)


class OpenMMSimulation(ParallelOEMolComputeCube):
    title = "OpenMM MD Simulation"
    version = "0.0.2"
    classification = [["Simulation", "OpenMM", "Minimization"],
    ["Simulation", "OpenMM", "Molecular Dynamics"]]
    tags = ['OpenMM', 'MDTraj', 'Parallel Cube']
    description = """
    Run an OpenMM molecular dynamics simulation. Default: 500K MD steps (1ns).

    Minimizes the prepared protein:ligand complex or restarts the simulation
    if there is an attached openmm.State on the molecule. Updates the attached
    parmed.Structure from the final state of the simulation. Attaches the
    openmm.State and output log file from the simulation energy reporter.

    Input:
    -------
    oechem.OEMol - Requires a streamed-in 'packed' molecule containing:
        - Generic Tags: { Structure : parmed.Structure (base64-encoded) }

    Restarts requires:
        - Generic Tags: { Structure : parmed.Structure (base64-encoded),
                          State : openmm.State (serialized) }

    Output:
    -------
    oechem.OEMol - Emits molecule with attachments:
        - SDData Tags: { Structure: str <parmed.Structure> }
        - Generic Tags: { Structure : parmed.Structure (base64-encoded),
                          State : openmm.State (serialized),
                          Log : txt-file (serialized) }
    tarxz - Generated tarball (LZMA compressed) containing simulation data:
        - Energy log, State.XML, Trajectory (default: NetCDF), and a PDB.
    """
    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 28800}, # Default 8 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)")

    steps = parameter.IntegerParameter(
        'steps',
        default=500000,
        help_text="Number of MD steps (500K = 1ns)")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=10000,
        help_text="Step interval for reporting data.")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The nonbonded cutoff in angstroms.
        This is ignored if nonbondedMethod is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    trajectory_filetype = parameter.StringParameter(
        'trajectory_filetype',
        default='NetCDF',
        choices=['NetCDF', 'DCD', 'HDF5'],
        help_text="NetCDF, DCD, HDF5. Filetype to write trajectory files")

    trajectory_selection = parameter.StringParameter(
        'trajectory_selection',
        default=None,
        choices=[None, 'protein or resname LIG', 'protein', 'resname LIG'],
        help_text='atoms subset to write in trajectory')

    trajectory_interval = parameter.IntegerParameter(
        'trajectory_interval',
        default=1000,
        help_text="Step interval for trajetory snapshots.")

    outfname = parameter.StringParameter(
        'outfname',
        default='md',
        help_text='Filename suffix for output simulation files. Formatted: <title>-<outfname>')

    tarxz = parameter.BooleanParameter(
        'tarxz',
        default=True,
        description='Create a tar.xz file of the attached data')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity.')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['convert'] = False
        self.opt['Logger'] = self.log
        conv_rule = [self.opt['trajectory_selection'] != None,
                     self.opt['trajectory_filetype'] != 'NetCDF']

        if any(conv_rule):
            self.opt['convert'] = True

    def process(self, mol, port):
        try:
            # Check for generic data.
            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                self.opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            # Generate Simulation from Structure
            simulation = simtools.genSimFromStruct(gd['Structure'], **self.opt)

            # Check if mol has State data attached
            if 'State' in gd.keys():
                self.log.info('%s RESTARTING from saved State' % gd['IDTag'])
                simulation.context.setState(gd['State'])
            else:
                self.log.info('%s MINIMIZING System' % gd['IDTag'])
                minene, simulation = simtools.minimizeSimulation(simulation, **self.opt)
                oechem.OESetSDData(mol, 'Minimized Energy', str(minene))

            for rep in simtools.getReporters(**self.opt):
                simulation.reporters.append(rep)

            self.log.info('{} running {steps} MD steps at {temperature}K'.format(
                gd['IDTag'], **self.opt))
            simulation.step(self.args.steps)

            if self.opt['convert']:
                self.log.info(
                    'Converting trajectories to: {trajectory_filetype}'.format(**self.opt))
                simtools.mdTrajConvert(simulation, outfname=self.opt['outfname'],
                                       trajectory_selection=self.opt['trajectory_selection'],
                                       trajectory_filetype=self.opt['trajectory_filetype'])

            packedmol = utils.PackageOEMol.pack(mol, simulation)
            packedmol.SetData(oechem.OEGetTag(
                'outfname'), self.opt['outfname'])

            # Create a tar.xz archive of the generic data and trajectories
            if self.opt['tarxz']:
                utils.PackageOEMol.dump(
                    packedmol, outfname=self.opt['outfname'], tarxz=self.opt['tarxz'])
            self.success.emit(packedmol)

        except Exception as e:
                # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)
