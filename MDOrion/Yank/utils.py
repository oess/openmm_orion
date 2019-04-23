# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


from yank.experiment import ExperimentBuilder

from MDOrion.Yank.yank_templates import (max_cube_running_time,
                                         yank_binding_template,
                                         yank_solvation_template)

from yank.analyze import ExperimentAnalyzer

from simtk.openmm import unit

from simtk import openmm

from datetime import timedelta

import subprocess

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

import yaml

import time

import fcntl

import os

from openeye import (oedepict,
                     oechem)

from tempfile import TemporaryDirectory


def yank_solvation_initialize(sim):
    def wrapper(*args):

        opt = args[0]

        yank_template = yank_solvation_template.format(
            verbose='yes' if opt['verbose'] else 'no',
            minimize='yes' if opt['minimize'] else 'no',
            output_directory=opt['out_directory'],
            timestep=4.0 if opt['hmr'] else 2.0,
            nsteps_per_iteration=opt['nsteps_per_iteration'],
            number_iterations=opt['new_iterations'],
            temperature=opt['temperature'],
            pressure=opt['pressure'],
            resume_sim='yes' if opt['resume_sim'] else 'no',
            resume_setup='yes' if opt['resume_setup'] else 'no',
            hydrogen_mass=4.0 if opt['hmr'] else 1.0,
            alchemical_pme_treatment=opt['alchemical_pme_treatment'],
            checkpoint_interval=opt['checkpoint_interval'],
            solvated_pdb_fn=opt['solvated_structure_fn'],
            solvated_xml_fn=opt['solvated_omm_sys_serialized_fn'],
            solute_pdb_fn=opt['solute_structure_fn'],
            solute_xml_fn=opt['solute_omm_sys_serialized_fn'],
            solvent_dsl=opt['solvent_str_names'])

        opt['yank_template'] = yank_template

        # Print Yank Template
        opt['Logger'].info(opt['yank_template'])

        if 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(',')

            opt['Logger'].info("OE LOCAL FLOE CLUSTER OPTION IN USE")

            if 'OE_MAX' in os.environ:
                opt['OE_MAX'] = int(os.environ["OE_MAX"])
            else:
                opt['OE_MAX'] = 1

            opt['Logger'].info("OE MAX = {}".format(opt['OE_MAX']))

            while True:

                for gpu_id in gpus_available_indexes:

                    for p in range(0, opt['OE_MAX']):

                        fn = str(gpu_id) + '_' + str(p) + '.txt'

                        try:
                            with open(fn, 'a') as file:
                                fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                                file.write(
                                    "YANK - name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(opt['system_title'],
                                                                                                     opt['system_id'],
                                                                                                     gpus_available_indexes,
                                                                                                     str(gpu_id)))
                                opt['gpu_id'] = str(gpu_id)

                                openmm.Platform.getPlatformByName('CUDA').setPropertyDefaultValue('DeviceIndex',
                                                                                                  str(gpu_id))

                                sim(opt)

                                time.sleep(5.0)
                                fcntl.flock(file, fcntl.LOCK_UN)
                                return

                        except BlockingIOError:
                            time.sleep(0.1)

                        except Exception as e:  # If the simulation fails for other reasons
                            try:
                                time.sleep(5.0)
                                fcntl.flock(file, fcntl.LOCK_UN)
                            except Exception as e:
                                pass
                            raise ValueError("{} Simulation Failed".format(str(e)))
        else:
            sim(*args)

    return wrapper


def yank_binding_initialize(sim):
    def wrapper(*args):

        opt = args[0]

        yank_template = yank_binding_template.format(
            verbose='yes' if opt['verbose'] else 'no',
            minimize='yes' if opt['minimize'] else 'no',
            output_directory=opt['out_directory'],
            timestep=opt['timestep'],
            nsteps_per_iteration=opt['nsteps_per_iteration'],
            number_iterations=opt['new_iterations'],
            temperature=opt['temperature'],
            pressure=opt['pressure'],
            resume_sim='yes' if opt['resume_sim'] else 'no',
            resume_setup='yes' if opt['resume_setup'] else 'no',
            hydrogen_mass=4.0 if opt['hmr'] else 1.0,
            alchemical_pme_treatment=opt['alchemical_pme_treatment'],
            checkpoint_interval=opt['checkpoint_interval'],
            complex_pdb_fn=opt['solvated_complex_structure_fn'],
            complex_xml_fn=opt['solvated_complex_omm_serialized_fn'],
            solvent_pdb_fn=opt['solvated_ligand_structure_fn'],
            solvent_xml_fn=opt['solvated_ligand_omm_serialized_fn'],
            ligand_resname=opt['lig_res_name'],
            solvent_dsl=opt['solvent_str_names'],
            sampler=opt['sampler'],
            restraints=opt['restraints'],
            protocol=opt['protocol'])

        if opt['user_yank_yaml_file'] is not None:

            opt['Logger'].info("Yank User template file in use ")

            # Loading the yank_template as yaml
            yank_yaml_template = yaml.load(yank_template)

            # Loading the user Yaml file
            try:
                if in_orion():
                    if not isinstance(opt['user_yank_yaml_file'], dict):
                        fn = opt['user_yank_yaml_file']
                    else:
                        file_id = opt['user_yank_yaml_file']['file']
                        session = OrionSession()
                        resource = session.get_resource(File, file_id)
                        fn = os.path.join(opt['output_directory'], "user_yank_orion.yaml")
                        resource.download_to_file(fn)
                else:
                    fn = opt['user_yank_yaml_file']

                with open(fn, 'r') as yaml_file:
                    yank_yaml_user = yaml.load(yaml_file)
            except:
                raise IOError("It was not possible to load the provided yaml file {}".format(opt['user_yank_yaml_file']))

            yank_yaml_merge = dict(yank_yaml_template)

            # Cleaning
            del yank_yaml_merge['harmonic']
            del yank_yaml_merge['boresch']

            if 'experiments' in yank_yaml_user:
                if isinstance(yank_yaml_user['experiments'], list):

                    for exp in yank_yaml_user['experiments']:
                        yank_yaml_merge[exp] = yank_yaml_user[exp]

                    yank_yaml_merge['experiments'] = yank_yaml_user['experiments']

                else:
                    yank_yaml_merge['experiments'] = ["exp"]

                    yank_yaml_merge['exp'] = yank_yaml_user['experiments']

                if len(yank_yaml_merge['experiments']) > 1:
                    raise ValueError("Currently just one Yank experiment is supported. Provided: {}".format(yank_yaml_merge['experiments']))

                for exp in yank_yaml_merge['experiments']:

                    # System section
                    yank_yaml_merge[exp]['system'] = 'system'

                    # Protocol Section
                    if 'protocol' in yank_yaml_merge[exp]:
                        yank_yaml_merge['protocols'][yank_yaml_merge[exp]['protocol']] = yank_yaml_user['protocols'][yank_yaml_merge[exp]['protocol']]

                    else:  # If the protocol section is not present use the selected one
                        yank_yaml_merge[exp]['protocol'] = yank_yaml_template[yank_yaml_template['experiments'][0]]['protocol']

                    if 'sampler' in yank_yaml_merge[exp]:
                        yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']] = yank_yaml_user['samplers'][yank_yaml_merge[exp]['sampler']]
                        if 'mcmc_moves' in yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]:
                            yank_yaml_merge['mcmc_moves'][yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]['mcmc_moves']] = yank_yaml_user['mcmc_moves'][yank_yaml_user['samplers'][yank_yaml_merge[exp]['sampler']]['mcmc_moves']]
                        else:  # mcmc move is set to the default one
                            yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]['mcmc_moves'] = yank_yaml_merge['samplers'][yank_yaml_template[yank_yaml_template['experiments'][0]]['sampler']]['mcmc_moves']
                    else:  # The sampler is set to the default one
                        yank_yaml_merge[exp]['sampler'] = yank_yaml_template[yank_yaml_template['experiments'][0]]['sampler']

                    # Sampler section number of iterations
                    yank_yaml_merge['samplers'][yank_yaml_merge[exp]['sampler']]['number_of_iterations'] = opt['new_iterations']

                yank_template = yank_yaml_merge
            else:
                raise ValueError("Missing Experiments section in the provided Yank Yaml file")

        opt['yank_template'] = yank_template

        # Print Yank Template
        if isinstance(opt['yank_template'], dict):
            opt['Logger'].info(yaml.dump(opt['yank_template']))
        else:
            opt['Logger'].info(opt['yank_template'])

        if 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(',')

            opt['Logger'].info("OE LOCAL FLOE CLUSTER OPTION IN USE")

            if 'OE_MAX' in os.environ:
                opt['OE_MAX'] = int(os.environ["OE_MAX"])
            else:
                opt['OE_MAX'] = 1

            opt['Logger'].info("OE MAX = {}".format(opt['OE_MAX']))

            while True:

                for gpu_id in gpus_available_indexes:

                    for p in range(0, opt['OE_MAX']):

                        fn = str(gpu_id) + '_' + str(p) + '.txt'

                        try:
                            with open(fn, 'a') as file:
                                fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                                file.write(
                                    "YANK - name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(opt['system_title'],
                                                                                                     opt['system_id'],
                                                                                                     gpus_available_indexes,
                                                                                                     str(gpu_id)))
                                opt['gpu_id'] = str(gpu_id)

                                openmm.Platform.getPlatformByName('CUDA').setPropertyDefaultValue('DeviceIndex',
                                                                                                  str(gpu_id))

                                sim(opt)

                                time.sleep(5.0)
                                fcntl.flock(file, fcntl.LOCK_UN)
                                return
                        except BlockingIOError:
                            time.sleep(0.1)

                        except Exception as e:  # If the simulation fails for other reasons
                            try:
                                time.sleep(5.0)
                                fcntl.flock(file, fcntl.LOCK_UN)
                            except Exception as e:
                                pass
                            raise ValueError("{} Simulation Failed".format(str(e)))
        else:
            sim(*args)

    return wrapper


@yank_solvation_initialize
def run_yank_solvation(opt):

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return


@yank_binding_initialize
def run_yank_binding(opt):

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return


def run_yank_analysis(opt):

    exp_dir = os.path.join(opt['out_directory'], "experiments")

    experiment_to_analyze = ExperimentAnalyzer(exp_dir)
    analysis = experiment_to_analyze.auto_analyze()

    # Calculate free energy and its error
    DeltaG = analysis['free_energy']['free_energy_diff_unit'].in_units_of(
        unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
    dDeltaG = analysis['free_energy']['free_energy_diff_error_unit'].in_units_of(
        unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole

    opt_1 = '--store={}'.format(exp_dir)

    result_fn = os.path.join(opt['out_directory'], 'results.html')
    opt_2 = '--output={}'.format(result_fn)

    opt_3 = '--format=html'

    report_html_str = None

    try:
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        subprocess.check_call(['yank', 'analyze', 'report', opt_1, opt_2, opt_3])

        result_fn = os.path.join(opt['out_directory'], "results.html")

        with open(result_fn, 'r') as f:
            report_html_str = f.read()

    except subprocess.SubprocessError:
        opt['Logger'].warn("The result html file has not been generated")

    return DeltaG, dDeltaG, report_html_str


def calculate_iteration_time(output_directory, num_iterations):

    UNITS = {"s": "seconds", "m": "minutes", "h": "hours", "d": "days", "w": "weeks"}

    def convert_to_seconds(s):
        count = int(float(s[:-1]))
        user_unit = UNITS[s[-1]]
        td = timedelta(**{user_unit: count})
        return td.seconds + 60 * 60 * 24 * td.days

    log_fn = os.path.join(output_directory, "experiments/experiments.log")

    f = open(log_fn, 'r')

    log = f.read()

    times = []
    for line in log.splitlines():
        if "Iteration took" in line:
            times.append(line)

    clean_times = [i.split('-')[-1].split(" ")[-1][:-1] for i in times]

    # times in seconds
    tms_list = []

    for tm in clean_times:
        tms = convert_to_seconds(tm)
        tms_list.append(tms)

    # average time per iterations in hrs
    avg_time_per_iteration_phase_1 = (sum(tms_list[:num_iterations]) / num_iterations) / 3600.0
    avg_time_per_iteration_phase_2 = (sum(tms_list[num_iterations:]) / num_iterations) / 3600.0

    avg_time_per_iteration = avg_time_per_iteration_phase_1 + avg_time_per_iteration_phase_2

    iterations_per_cube = int(max_cube_running_time / avg_time_per_iteration)

    return iterations_per_cube


def ligand_to_svg(ligand, ligand_name):

    with TemporaryDirectory() as output_directory:

        lig_copy = oechem.OEMol(ligand)

        if len(ligand_name) < 15:
            lig_copy.SetTitle(ligand_name)
        else:
            lig_copy.SetTitle(ligand_name[0:13] + '...')

        img_fn = os.path.join(output_directory, "img.svg")

        oedepict.OEPrepareDepiction(lig_copy)

        width, height = 150, 150
        opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
        opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)
        disp = oedepict.OE2DMolDisplay(lig_copy, opts)

        oedepict.OERenderMolecule(img_fn, disp)

        svg_lines = ""
        marker = False
        with open(img_fn, 'r') as file:
            for line in file:
                if marker:
                    svg_lines += line

                if line.startswith("<svg"):
                    marker = True
                    svg_lines += line
                    svg_lines += """<title>{}</title>\n""".format(ligand_name)

                if line.startswith("</svg>"):
                    marker = False

    return svg_lines
