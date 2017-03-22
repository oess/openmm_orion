import io, os, time, traceback, base64, shutil
from openeye import oechem
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app
import netCDF4 as netcdf
from tempfile import TemporaryDirectory

from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort, ParallelOEMolComputeCube
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
import YankCubes.utils as utils
from YankCubes.utils import molecule_is_charged, download_dataset_to_file, get_data_filename

import yank
import mdtraj
from yank.yamlbuild import *
import textwrap
import subprocess

################################################################################
# Hydration free energy calculations
################################################################################

hydration_yaml_template = """\
---
options:
  minimize: no
  timestep: %(timestep)f*femtoseconds
  nsteps_per_iteration: %(nsteps_per_iteration)d
  number_of_iterations: %(number_of_iterations)d
  temperature: %(temperature)f*kelvin
  pressure: %(pressure)f*atmosphere
  anisotropic_dispersion_correction: yes
  anisotropic_dispersion_cutoff: 9*angstroms
  output_dir: %(output_directory)s
  verbose: %(verbose)s

molecules:
  input_molecule:
    filepath: %(output_directory)s/input.mol2
    antechamber:
      charge_method: null

solvents:
  tip3p:
    nonbonded_method: PME
    nonbonded_cutoff: 9*angstroms
    clearance: 8*angstroms
    ewald_error_tolerance: 1.0e-4
  gbsa:
    nonbonded_method: NoCutoff
    implicit_solvent: OBC2
  vacuum:
    nonbonded_method: NoCutoff

systems:
  hydration-tip3p:
    solute: input_molecule
    solvent1: tip3p
    solvent2: vacuum
    leap:
      parameters: [leaprc.gaff, leaprc.protein.ff14SB, leaprc.water.tip3p]
  hydration-gbsa:
    solute: input_molecule
    solvent1: gbsa
    solvent2: vacuum
    leap:
      parameters: [leaprc.gaff, leaprc.protein.ff14SB, leaprc.water.tip3p]

protocols:
  protocol-tip3p:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.00]
  protocol-gbsa:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.00]
        lambda_sterics:        [1.00, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.00]
        lambda_sterics:        [1.00, 0.00]

experiments:
  system: hydration-%(solvent)s
  protocol: protocol-%(solvent)s
"""

class YankHydrationCube(ParallelOEMolComputeCube):
    title = "YankHydrationCube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to TIP3P solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = ["Alchemical free energy calculations"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')
    failure = CustomMoleculeOutputPort('failure')

    # These can override YAML parameters
    nsteps_per_iteration = parameter.IntegerParameter('nsteps_per_iteration', default=500,
                                     help_text="Number of steps per iteration")

    timestep = parameter.DecimalParameter('timestep', default=2.0,
                                     help_text="Timestep (fs)")

    simulation_time = parameter.DecimalParameter('simulation_time', default=0.100,
                                     help_text="Simulation time (ns/replica)")

    temperature = parameter.DecimalParameter('temperature', default=300.0,
                                     help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter('pressure', default=1.0,
                                 help_text="Pressure (atm)")

    solvent = parameter.StringParameter('solvent', default='gbsa',
                                 choices=['gbsa', 'tip3p'],
                                 help_text="Solvent choice: one of ['gbsa', 'tip3p']")

    verbose = parameter.BooleanParameter('verbose', default=False,
                                     help_text="Print verbose YANK logging output")

    def construct_yaml(self, **kwargs):
        # Make substitutions to YAML here.
        # TODO: Can we override YAML parameters without having to do string substitutions?
        options = {
            'timestep' : self.args.timestep,
            'nsteps_per_iteration' : self.args.nsteps_per_iteration,
            'number_of_iterations' : int(np.ceil(self.args.simulation_time * unit.nanoseconds / (self.args.nsteps_per_iteration * self.args.timestep * unit.femtoseconds))),
            'temperature' : self.args.temperature,
            'pressure' : self.args.pressure,
            'solvent' : self.args.solvent,
            'verbose' : 'yes' if self.args.verbose else 'no',
        }

        for parameter in kwargs.keys():
            options[parameter] = kwargs[parameter]

        return hydration_yaml_template % options

    def begin(self):
        # TODO: Is there another idiom to use to check valid input?
        if self.args.solvent not in ['gbsa', 'tip3p']:
            raise Exception("solvent must be one of ['gbsa', 'tip3p']")

        # Compute kT
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA # Boltzmann constant
        self.kT = kB * (self.args.temperature * unit.kelvin)

    def process(self, mol, port):
        kT_in_kcal_per_mole = self.kT.value_in_unit(unit.kilocalories_per_mole)

        # Retrieve data about which molecule we are processing
        title = mol.GetTitle()

        with TemporaryDirectory() as output_directory:
            try:
                # Print out which molecule we are processing
                self.log.info('Processing {} in directory {}.'.format(title, output_directory))

                # Check that molecule is charged.
                if not molecule_is_charged(mol):
                    raise Exception('Molecule %s has no charges; input molecules must be charged.' % mol.GetTitle())

                # Write the specified molecule out to a mol2 file without changing its name.
                mol2_filename = os.path.join(output_directory, 'input.mol2')
                ofs = oechem.oemolostream(mol2_filename)
                oechem.OEWriteMol2File(ofs, mol)

                # Undo oechem fuckery with naming mol2 substructures `<0>`
                from YankCubes.utils import unfuck_oechem_mol2_file
                unfuck_oechem_mol2_file(mol2_filename)

                # Run YANK on the specified molecule.
                from yank.yamlbuild import YamlBuilder
                yaml = self.construct_yaml(output_directory=output_directory)
                yaml_builder = YamlBuilder(yaml)
                yaml_builder.build_experiments()
                self.log.info('Ran Yank experiments for molecule {}.'.format(title))

                # Analyze the hydration free energy.
                from yank.analyze import estimate_free_energies
                (Deltaf_ij_solvent, dDeltaf_ij_solvent) = estimate_free_energies(netcdf.Dataset(output_directory + '/experiments/solvent1.nc', 'r'))
                (Deltaf_ij_vacuum,  dDeltaf_ij_vacuum)  = estimate_free_energies(netcdf.Dataset(output_directory + '/experiments/solvent2.nc', 'r'))
                DeltaG_hydration = Deltaf_ij_vacuum[0,-1] - Deltaf_ij_solvent[0,-1]
                dDeltaG_hydration = np.sqrt(Deltaf_ij_vacuum[0,-1]**2 + Deltaf_ij_solvent[0,-1]**2)

                # Add result to original molecule
                oechem.OESetSDData(mol, 'DeltaG_yank_hydration', str(DeltaG_hydration * kT_in_kcal_per_mole))
                oechem.OESetSDData(mol, 'dDeltaG_yank_hydration', str(dDeltaG_hydration * kT_in_kcal_per_mole))
                self.log.info('Analyzed and stored hydration free energy for molecule {}.'.format(title))

                # Emit molecule to success port.
                self.success.emit(mol)

            except Exception as e:
                self.log.info('Exception encountered when processing molecule {}.'.format(title))
                # Attach error message to the molecule that failed
                # TODO: If there is an error in the leap setup log,
                # we should capture that and attach it to the failed molecule.
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed molecule
                self.failure.emit(mol)

################################################################################
# Binding free energy calculations
################################################################################

binding_yaml_template = """\
---
options:
  minimize: %(minimize)s
  timestep: %(timestep)f*femtoseconds
  nsteps_per_iteration: %(nsteps_per_iteration)d
  number_of_iterations: %(number_of_iterations)d
  temperature: %(temperature)f*kelvin
  pressure: %(pressure)f*atmosphere
  output_dir: %(output_directory)s
  verbose: %(verbose)s
  randomize_ligand: %(randomize_ligand)s

molecules:
  receptor:
    filepath: %(output_directory)s/receptor.pdb
    strip_protons: yes
  ligand:
    filepath: %(output_directory)s/input.mol2
    antechamber:
      charge_method: null

solvents:
  pme:
    nonbonded_method: PME
    nonbonded_cutoff: 9*angstroms
    ewald_error_tolerance: 1.0e-4
    clearance: 9*angstroms
    positive_ion: Na+
    negative_ion: Cl-
  rf:
    nonbonded_method: CutoffPeriodic
    nonbonded_cutoff: 9*angstroms
    clearance: 9*angstroms
    positive_ion: Na+
    negative_ion: Cl-
  gbsa:
    nonbonded_method: NoCutoff
    implicit_solvent: OBC2

systems:
  binding-pme:
    receptor: receptor
    ligand: ligand
    solvent: pme
    leap:
      parameters: [leaprc.protein.ff14SB, leaprc.gaff2, leaprc.water.tip3p]
  binding-rf:
    receptor: receptor
    ligand: ligand
    solvent: rf
    leap:
      parameters: [leaprc.protein.ff14SB, leaprc.gaff2, leaprc.water.tip3p]
  binding-gbsa:
    receptor: receptor
    ligand: ligand
    solvent: gbsa
    leap:
      parameters: [leaprc.protein.ff14SB, leaprc.gaff2, leaprc.water.tip3p]

protocols:
  protocol:
    complex:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]
        lambda_restraints:     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
    solvent:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]

experiments:
  system: binding-%(solvent)s
  protocol: protocol
  restraint:
    type: Harmonic
"""

class YankBindingCube(ParallelOEMolComputeCube):
    title = "YankBindingCube"
    description = """
    Compute thebinding free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the binding
    free energy of one or more small molecules using harmonic restraints.

    See http://getyank.org for more information about YANK.
    """
    classification = ["Alchemical free energy calculations"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')
    failure = CustomMoleculeOutputPort('failure')

    # Receptor specification
    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Receptor structure file')

    # These can override YAML parameters
    nsteps_per_iteration = parameter.IntegerParameter('nsteps_per_iteration', default=500,
                                     help_text="Number of steps per iteration")

    timestep = parameter.DecimalParameter('timestep', default=2.0,
                                     help_text="Timestep (fs)")

    simulation_time = parameter.DecimalParameter('simulation_time', default=0.100,
                                     help_text="Simulation time (ns/replica)")

    temperature = parameter.DecimalParameter('temperature', default=300.0,
                                     help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter('pressure', default=1.0,
                                 help_text="Pressure (atm)")

    solvent = parameter.StringParameter('solvent', default='gbsa',
                                 choices=['gbsa', 'pme', 'rf'],
                                 help_text="Solvent choice ['gbsa', 'pme', 'rf']")

    minimize = parameter.BooleanParameter('minimize', default=True,
                                     help_text="Minimize initial structures for stability")

    randomize_ligand = parameter.BooleanParameter('randomize_ligand', default=False,
                                     help_text="Randomize initial ligand position (implicit only)")

    verbose = parameter.BooleanParameter('verbose', default=False,
                                     help_text="Print verbose YANK logging output")

    def construct_yaml(self, **kwargs):
        # Make substitutions to YAML here.
        # TODO: Can we override YAML parameters without having to do string substitutions?
        options = {
            'timestep' : self.args.timestep,
            'nsteps_per_iteration' : self.args.nsteps_per_iteration,
            'number_of_iterations' : int(np.ceil(self.args.simulation_time * unit.nanoseconds / (self.args.nsteps_per_iteration * self.args.timestep * unit.femtoseconds))),
            'temperature' : self.args.temperature,
            'pressure' : self.args.pressure,
            'solvent' : self.args.solvent,
            'minimize' : 'yes' if self.args.minimize else 'no',
            'verbose' : 'yes' if self.args.verbose else 'no',
            'randomize_ligand' : 'yes' if self.args.randomize_ligand else 'no',
        }

        for parameter in kwargs.keys():
            options[parameter] = kwargs[parameter]

        return binding_yaml_template % options

    def begin(self):
        # TODO: Is there another idiom to use to check valid input?
        if self.args.solvent not in ['gbsa', 'pme', 'rf']:
            raise Exception("solvent must be one of ['gbsa', 'pme', 'rf']")

        # Compute kT
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA # Boltzmann constant
        self.kT = kB * (self.args.temperature * unit.kelvin)

        # Load receptor
        self.receptor = oechem.OEMol()
        receptor_filename = download_dataset_to_file(self.args.receptor)
        with oechem.oemolistream(receptor_filename) as ifs:
            if not oechem.OEReadMolecule(ifs, self.receptor):
                raise RuntimeError("Error reading receptor")

    def process(self, mol, port):
        kT_in_kcal_per_mole = self.kT.value_in_unit(unit.kilocalories_per_mole)

        # Retrieve data about which molecule we are processing
        title = mol.GetTitle()

        with TemporaryDirectory() as output_directory:
            try:
                # Print out which molecule we are processing
                self.log.info('Processing {} in {}.'.format(title, output_directory))

                # Check that molecule is charged.
                if not molecule_is_charged(mol):
                    raise Exception('Molecule %s has no charges; input molecules must be charged.' % mol.GetTitle())


                # Write the receptor.
                pdbfilename = os.path.join(output_directory, 'receptor.pdb')
                with oechem.oemolostream(pdbfilename) as ofs:
                    res = oechem.OEWriteConstMolecule(ofs, self.receptor)
                    if res != oechem.OEWriteMolReturnCode_Success:
                        raise RuntimeError("Error writing receptor: {}".format(res))

                # Write the specified molecule out to a mol2 file without changing its name.
                mol2_filename = os.path.join(output_directory, 'input.mol2')
                ofs = oechem.oemolostream(mol2_filename)
                oechem.OEWriteMol2File(ofs, mol)

                # Undo oechem fuckery with naming mol2 substructures `<0>`
                from YankCubes.utils import unfuck_oechem_mol2_file
                unfuck_oechem_mol2_file(mol2_filename)

                # Run YANK on the specified molecule.
                from yank.yamlbuild import YamlBuilder
                yaml = self.construct_yaml(output_directory=output_directory)
                yaml_builder = YamlBuilder(yaml)
                yaml_builder.build_experiments()
                self.log.info('Ran Yank experiments for molecule {}.'.format(title))

                # Analyze the binding free energy
                # TODO: Use yank.analyze API for this
                from YankCubes.analysis import analyze
                store_directory = os.path.join(output_directory, 'experiments')
                [DeltaG_binding, dDeltaG_binding] = analyze(store_directory)

                """
                # Extract trajectory (DEBUG)
                from yank.analyze import extract_trajectory
                trajectory_filename = 'trajectory.pdb'
                store_filename = os.path.join(store_directory, 'complex.pdb')
                extract_trajectory(trajectory_filename, store_filename, state_index=0, keep_solvent=False,
                       discard_equilibration=True, image_molecules=True)
                ifs = oechem.oemolistream(trajectory_filename)
                ifs.SetConfTest(oechem.OEAbsCanonicalConfTest()) # load multi-conformer molecule
                mol = oechem.OEMol()
                for mol in ifs.GetOEMols():
                    print (mol.GetTitle(), "has", mol.NumConfs(), "conformers")
                ifs.close()
                os.remove(trajectory_filename)
                """

                # Attach binding free energy estimates to molecule
                oechem.OESetSDData(mol, 'DeltaG_yank_binding', str(DeltaG_binding * kT_in_kcal_per_mole))
                oechem.OESetSDData(mol, 'dDeltaG_yank_binding', str(dDeltaG_binding * kT_in_kcal_per_mole))
                self.log.info('Analyzed and stored binding free energy for molecule {}.'.format(title))

                # Emit molecule to success port.
                self.success.emit(mol)

            except Exception as e:
                self.log.info('Exception encountered when processing molecule {}.'.format(title))
                # Attach error message to the molecule that failed
                # TODO: If there is an error in the leap setup log,
                # we should capture that and attach it to the failed molecule.
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed molecule
                self.failure.emit(mol)
