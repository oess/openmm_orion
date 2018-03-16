import traceback
from openeye import oechem
from tempfile import TemporaryDirectory
from floe.api import (parameter, ParallelMixin, ParallelOEMolComputeCube, OEMolComputeCube, MoleculeInputPort,
                      BatchMoleculeOutputPort, BatchMoleculeInputPort)
from floe.api.orion import in_orion

from cuberecord import OERecordComputeCube, OEField
from cuberecord.constants import DEFAULT_MOL_NAME
from datarecord import Types, Meta, ColumnMeta

from yank.experiment import ExperimentBuilder
from oeommtools import utils as oeommutils
from oeommtools import data_utils
from simtk.openmm import app, unit, XmlSerializer, openmm
import os
from YankCubes import utils as yankutils
from YankCubes.yank_templates import yank_solvation_template, yank_binding_template
import itertools
import OpenMMCubes.utils as omm_utils

import tarfile
from big_storage import LargeFileDataType

from mdtraj.core.residue_names import _WATER_RESIDUES as water_names


class YankSolvationFECube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "YankSolvationFECube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to the selected solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = ["Alchemical free energy calculations"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda8"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    minimize = parameter.BooleanParameter(
        'minimize',
        default=True,
        help_text="Minimize input system")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=1000,
        help_text="Number of iterations")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=500,
        help_text="Number of steps per iteration")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10.0,
        help_text="The non-bonded cutoff in angstroms")

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text="Print verbose YANK logging output")

    rerun = parameter.BooleanParameter(
        'rerun',
        default=False,
        help_text="Start Yank Restart procedure")

    analyze = parameter.BooleanParameter(
        'analyze',
        default=False,
        help_text="Start Yank Analysis on the collected results")

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='Hydrogen Mass Repartitioning')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            field_system = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol,
                                   meta=ColumnMeta().set_option(Meta.Hints.Chem.PrimaryMol))

            if not record.has_value(field_system):
                self.log.warn("Missing molecule '{}' field".format(field_system.get_name()))
                self.failure.emit(record)
                return

            solvated_system = record.get_value(field_system)

            field_parmed = OEField("Parmed", omm_utils.ParmedData)

            if not record.has_value(field_parmed):
                self.log.warn("Missing molecule '{}' field".format(field_parmed.get_name()))
                self.failure.emit(record)

            parmed_structure = record.get_value(field_parmed)

            # Split the complex in components
            if not opt['rerun']:
                protein, solute, water, excipients = oeommutils.split(solvated_system, ligand_res_name='LIG')
            else:
                solute = solvated_system

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solute) if dp.GetTag() in
                        ["temperature", "pressure"]}

            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solute.GetTitle(), new_args))
                opt.update(new_args)

            # Extract the MD data
            mdData = omm_utils.MDData(parmed_structure)
            solvated_structure = mdData.structure

            # Extract the ligand parmed structure
            solute_structure = solvated_structure.split()[0][0]
            solute_structure.box = None

            solvent_res_names = set()
            for res in solvated_structure.residues:
                solvent_res_names.add(res.name)
            solvent_res_names.remove(solute_structure.residues[0].name)

            solvent_str_names = ' '.join(solvent_res_names)

            # Testing
            solute_key = ''

            if not opt['rerun']:

                # Set the ligand title
                solute.SetTitle(solvated_system.GetTitle())

                # Create the solvated and vacuum system
                solvated_omm_sys = solvated_structure.createSystem(nonbondedMethod=app.PME,
                                                                   nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                   constraints=app.HBonds,
                                                                   removeCMMotion=False)

                solute_omm_sys = solute_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                               constraints=app.HBonds,
                                                               removeCMMotion=False)

                # This is a note from:
                # https://github.com/MobleyLab/SMIRNOFF_paper_code/blob/e5012c8fdc4570ca0ec750f7ab81dd7102e813b9/scripts/create_input_files.py#L114
                # Fix switching function.
                for force in solvated_omm_sys.getForces():
                    if isinstance(force, openmm.NonbondedForce):
                        force.setUseSwitchingFunction(True)
                        force.setSwitchingDistance((opt['nonbondedCutoff'] - 1.0) * unit.angstrom)

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                if opt['rerun']:
                    if in_orion():
                        lf_file = OEField("lf_field", LargeFileDataType)
                    else:
                        lf_file = OEField("lf_field", Types.String)

                    file_id = record.get_value(lf_file)
                    filename = yankutils.download(file_id)

                    with tarfile.open(filename) as tar:
                        tar.extractall(path=output_directory)
                        # os.remove(filename)

                    # Disable minimization if restart is enabled
                    opt['minimize'] = False

                solvated_structure_fn = os.path.join(output_directory, "solvated.pdb")
                solute_structure_fn = os.path.join(output_directory, "solute.pdb")

                solvated_omm_sys_serialized_fn = os.path.join(output_directory, "solvated.xml")
                solute_omm_sys_serialized_fn = os.path.join(output_directory, "solute.xml")

                if not opt['rerun']:
                    solvated_structure.save(solvated_structure_fn, overwrite=True)
                    solute_structure.save(solute_structure_fn, overwrite=True)

                    solvated_omm_sys_serialized = XmlSerializer.serialize(solvated_omm_sys)
                    solvated_f = open(solvated_omm_sys_serialized_fn, 'w')
                    solvated_f.write(solvated_omm_sys_serialized)
                    solvated_f.close()

                    solute_omm_sys_serialized = XmlSerializer.serialize(solute_omm_sys)
                    solute_f = open(solute_omm_sys_serialized_fn, 'w')
                    solute_f.write(solute_omm_sys_serialized)
                    solute_f.close()

                self.log.warn(yank_solvation_template.format(
                                                 verbose='yes' if opt['verbose'] else 'no',
                                                 minimize='yes' if opt['minimize'] else 'no',
                                                 output_directory=output_directory,
                                                 timestep=4.0 if opt['hmr'] else 2.0,
                                                 nsteps_per_iteration=opt['nsteps_per_iteration'],
                                                 number_iterations=opt['iterations'],
                                                 temperature=opt['temperature'],
                                                 pressure=opt['pressure'],
                                                 resume_sim='yes' if opt['rerun'] else 'no',
                                                 resume_setup='yes' if opt['rerun'] else 'no',
                                                 hydrogen_mass=4.0 if opt['hmr'] else 1.0,
                                                 solvated_pdb_fn=solvated_structure_fn,
                                                 solvated_xml_fn=solvated_omm_sys_serialized_fn,
                                                 solute_pdb_fn=solute_structure_fn,
                                                 solute_xml_fn=solute_omm_sys_serialized_fn,
                                                 solvent_dsl=solvent_str_names))
                # Build the Yank Experiment

                yaml_builder = ExperimentBuilder(yank_solvation_template.format(
                                                 verbose='yes' if opt['verbose'] else 'no',
                                                 minimize='yes' if opt['minimize'] else 'no',
                                                 output_directory=output_directory,
                                                 timestep=4.0 if opt['hmr'] else 2.0,
                                                 nsteps_per_iteration=opt['nsteps_per_iteration'],
                                                 number_iterations=opt['iterations'],
                                                 temperature=opt['temperature'],
                                                 pressure=opt['pressure'],
                                                 resume_sim='yes' if opt['rerun'] else 'no',
                                                 resume_setup='yes' if opt['rerun'] else 'no',
                                                 hydrogen_mass=4.0 if opt['hmr'] else 1.0,
                                                 solvated_pdb_fn=solvated_structure_fn,
                                                 solvated_xml_fn=solvated_omm_sys_serialized_fn,
                                                 solute_pdb_fn=solute_structure_fn,
                                                 solute_xml_fn=solute_omm_sys_serialized_fn,
                                                 solute=solute_key,
                                                 solvent_dsl=solvent_str_names))

                # Run Yank
                yaml_builder.run_experiments()

                # Tar the temp dir with its content:
                tar_fn = os.path.basename(output_directory) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                if in_orion():
                    lf_file = OEField("lf_field", LargeFileDataType)
                else:
                    lf_file = OEField("lf_field", Types.String)

                lf = yankutils.upload(tar_fn)

                record.set_value(lf_file, lf)

                if opt['analyze']:

                    exp_dir = os.path.join(output_directory, "experiments")

                    # Calculate solvation free energy, solvation Enthalpy and their errors
                    DeltaG_solvation, dDeltaG_solvation, DeltaH, dDeltaH = yankutils.analyze_directory(exp_dir)

                    # # Add result to the original molecule in kcal/mol
                    oechem.OESetSDData(solute, 'DG_yank_solv', str(DeltaG_solvation))
                    oechem.OESetSDData(solute, 'dG_yank_solv', str(dDeltaG_solvation))

            record.set_value(field_system, solute)

            # Emit the ligand
            self.success.emit(record)

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class SyncBindingFECube(OEMolComputeCube):
    version = "0.0.0"
    title = "SyncSolvationFECube"
    description = """
    This cube is used to synchronize the solvated ligands and the related
    solvated complexes 
    """
    classification = ["Synchronization Cube"]
    tags = [tag for lists in classification for tag in lists]

    solvated_ligand_in_port = MoleculeInputPort("solvated_ligand_in_port")

    # Define a molecule batch port to stream out the solvated ligand and complex
    solvated_lig_complex_out_port = BatchMoleculeOutputPort("solvated_lig_complex_out_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.solvated_ligand_list = []
        self.solvated_complex_list = []

    def process(self, solvated_system, port):

        try:
            if port == 'solvated_ligand_in_port':
                self.solvated_ligand_list.append(solvated_system)
            else:
                self.solvated_complex_list.append(solvated_system)

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            solvated_system.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(solvated_system)

        return

    def end(self):
        solvated_complex_lig_list = [(i, j) for i, j in
                                     itertools.product(self.solvated_ligand_list, self.solvated_complex_list)
                                     if i.GetData("IDTag") in j.GetData("IDTag")]

        for pair in solvated_complex_lig_list:
            print(pair[0].GetData("IDTag"), pair[1].GetData("IDTag"))
            self.solvated_lig_complex_out_port.emit([pair[0], pair[1]])

        return


class YankBindingFECube(ParallelOEMolComputeCube):
    version = "0.0.0"
    title = "YankSolvationFECube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to the selected solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = ["Alchemical free energy calculations"]
    tags = [tag for lists in classification for tag in lists]

    # The intake port is re-defined as batch port
    intake = BatchMoleculeInputPort("intake")

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    minimize = parameter.BooleanParameter(
        'minimize',
        default=False,
        help_text="Minimize input system")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=1000,
        help_text="Number of iterations")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=500,
        help_text="Number of steps per iteration")

    timestep = parameter.DecimalParameter(
        'timestep',
        default=2.0,
        help_text="Timestep (fs)")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10.0,
        help_text="The non-bonded cutoff in angstroms")

    restraints = parameter.StringParameter(
        'restraints',
        default='Harmonic',
        choices=['FlatBottom', 'Harmonic', 'Boresch'],
        help_text='Select the restraint types')

    ligand_resname = parameter.StringParameter(
        'ligand_resname',
        default='LIG',
        help_text='The decoupling ligand residue name')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        help_text="Print verbose YANK logging output")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, solvated_system, port):

        try:
            opt = dict(self.opt)

            # Extract the solvated ligand and the solvated complex
            solvated_ligand = solvated_system[0]
            solvated_complex = solvated_system[1]

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solvated_ligand) if dp.GetTag() in
                        ["temperature", "pressure"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solvated_ligand.GetTitle(), new_args))
                opt.update(new_args)

            # Extract the MD data
            mdData_ligand = data_utils.MDData(solvated_ligand)
            solvated_ligand_structure = mdData_ligand.structure

            mdData_complex = data_utils.MDData(solvated_complex)
            solvated_complex_structure = mdData_complex.structure

            # Create the solvated OpenMM systems
            solvated_complex_omm_sys = solvated_complex_structure.createSystem(nonbondedMethod=app.PME,
                                                                               nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                               constraints=app.HBonds,
                                                                               removeCMMotion=False)

            solvated_ligand_omm_sys = solvated_ligand_structure.createSystem(nonbondedMethod=app.PME,
                                                                             nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                             constraints=app.HBonds,
                                                                             removeCMMotion=False)

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                solvated_complex_structure_fn = os.path.join(output_directory, "complex.pdb")
                solvated_complex_structure.save(solvated_complex_structure_fn, overwrite=True)

                solvated_ligand_structure_fn = os.path.join(output_directory, "solvent.pdb")
                solvated_ligand_structure.save(solvated_ligand_structure_fn, overwrite=True)

                solvated_complex_omm_serialized = XmlSerializer.serialize(solvated_complex_omm_sys)
                solvated_complex_omm_serialized_fn = os.path.join(output_directory, "complex.xml")
                solvated_complex_f = open(solvated_complex_omm_serialized_fn, 'w')
                solvated_complex_f.write(solvated_complex_omm_serialized)
                solvated_complex_f.close()

                solvated_ligand_omm_serialized = XmlSerializer.serialize(solvated_ligand_omm_sys)
                solvated_ligand_omm_serialized_fn = os.path.join(output_directory, "solvent.xml")
                solvated_ligand_f = open(solvated_ligand_omm_serialized_fn, 'w')
                solvated_ligand_f.write(solvated_ligand_omm_serialized)
                solvated_ligand_f.close()

                # Build the Yank Experiment
                yaml_builder = ExperimentBuilder(yank_binding_template.format(
                    verbose='yes' if opt['verbose'] else 'no',
                    minimize='yes' if opt['minimize'] else 'no',
                    output_directory=output_directory,
                    timestep=opt['timestep'],
                    nsteps_per_iteration=opt['nsteps_per_iteration'],
                    number_iterations=opt['iterations'],
                    temperature=opt['temperature'],
                    pressure=opt['pressure'],
                    complex_pdb_fn=solvated_complex_structure_fn,
                    complex_xml_fn=solvated_complex_omm_serialized_fn,
                    solvent_pdb_fn=solvated_ligand_structure_fn,
                    solvent_xml_fn=solvated_ligand_omm_serialized_fn,
                    restraints=opt['restraints'],
                    ligand_resname=opt['ligand_resname']))

                # Run Yank
                yaml_builder.run_experiments()

                exp_dir = os.path.join(output_directory, "experiments")

                DeltaG_binding, dDeltaG_binding, DeltaH, dDeltaH = yankutils.analyze_directory(exp_dir)

                protein, ligand, water, excipients = oeommutils.split(solvated_ligand,
                                                                      ligand_res_name=opt['ligand_resname'])
                # Add result to the extracted ligand in kcal/mol
                oechem.OESetSDData(ligand, 'DG_yank_binding', str(DeltaG_binding))
                oechem.OESetSDData(ligand, 'dG_yank_binding', str(dDeltaG_binding))

            self.success.emit(ligand)

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            solvated_system[1].SetData('error', str(e))
            # Return failed mol
            self.failure.emit(solvated_system[1])

        return 