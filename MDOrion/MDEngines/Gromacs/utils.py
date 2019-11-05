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

import subprocess

from subprocess import STDOUT, PIPE, Popen, DEVNULL

import tempfile

import tarfile

import os

from MDOrion.MDEngines.Gromacs.standards import Gromacs

from orionclient.session import in_orion, APISession

from orionclient.types import File

from os import environ

from MDOrion.MDEngines.Gromacs.standards import Fields


def upload_file(filename, orion_ui_name='OrionFile', tag='Trajectory'):

    if in_orion():

        session = APISession

        file_upload = File.upload(session,
                                  orion_ui_name,
                                  filename)

        session.tag_resource(file_upload, tag)

        job_id = environ.get('ORION_JOB_ID')

        if job_id:
            session.tag_resource(file_upload, "Job {}".format(job_id))

        file_id = file_upload.id

    else:
        file_id = filename

    return file_id


def download_file(file_id, filename):

    if in_orion():

        session = APISession

        resource = session.get_resource(File, file_id)

        resource.download_to_file(filename)

        fn_local = filename

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The File has not been found: {}".format(fn_local))

    return fn_local


def gmx_steps(flag, filename):

    out = subprocess.check_output(["gmx",
                                   "dump",
                                   flag, filename])

    fn_dump_str = out.decode("utf-8")

    list_dump = fn_dump_str.split("\n")

    for l in list_dump:
        if 'step' in l:
            nsteps = int(l.split()[2])
            break

    return nsteps


def gmx_run(opt):

    cycle_id = opt['cycle_id']
    record = opt['record']

    prefix = opt['prefix_name']

    cwd = os.getcwd()

    if opt['cube_run_time'] > 12.0:
        raise ValueError("Max Cube Running time supported is 12 hrs. {} provided".format(opt['cube_run_time']))

    with tempfile.TemporaryDirectory() as outdir:

        # Write the .tpr file
        fn_tpr = os.path.join(outdir, prefix + "_" + Gromacs.gmx_tpr_prefix + '.tpr')

        f = open(fn_tpr, 'wb')
        f.write(opt['tpr'])
        f.close()

        restart_fn = os.path.join(outdir, prefix + "_" + Gromacs.gmx_restart_prefix + '_' + str(cycle_id) + '.cpt')

        if cycle_id == 0:

            opt['Logger'].info("....Initializing Iterations")

            os.chdir(outdir)

            # Run Gromacs
            if opt['verbose']:
                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', fn_tpr,
                                       '-deffnm', prefix + "_" + Gromacs.gmx_run_prefix,
                                       '-cpo', restart_fn,
                                       '-maxh', str(opt['cube_run_time'])
                                       ])
            else:
                p = Popen(['gmx',
                           'mdrun',
                           '-v',
                           '-s', fn_tpr,
                           '-deffnm', prefix + "_" + Gromacs.gmx_run_prefix,
                           '-cpo', restart_fn,
                           '-maxh', str(opt['cube_run_time'])
                           ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)
                p.communicate()

        else:

            opt['Logger'].info("....Restarting")

            # Download the Restart File from S3
            if not record.has_value(Fields.gmx_restart):
                raise ValueError("Missing the GMX restart Field")

            file_restart_id = record.get_value(Fields.gmx_restart)

            tar_restart_fn = prefix + "_" + Gromacs.gmx_restart_prefix + '_' + str(cycle_id) + '.tar'

            file_tar = download_file(file_restart_id, tar_restart_fn)

            with tarfile.open(file_tar) as tar:
                tar.extractall(path=outdir)

            restart_fn_prev = os.path.join(outdir, prefix + "_" + Gromacs.gmx_restart_prefix + '_' + str(cycle_id - 1) + '.cpt')

            os.chdir(outdir)

            # Run Gromacs
            if opt['verbose']:

                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', fn_tpr,
                                       '-deffnm', prefix + "_" + Gromacs.gmx_run_prefix,
                                       '-noappend',
                                       '-cpi', restart_fn_prev,
                                       '-cpo', restart_fn,
                                       '-maxh',  str(opt['cube_run_time'])
                                       ])
            else:
                p = Popen(['gmx',
                           'mdrun',
                           '-v',
                           '-s', fn_tpr,
                           '-deffnm', prefix + "_" + Gromacs.gmx_run_prefix,
                           '-noappend',
                           '-cpi', restart_fn_prev,
                           '-cpo', restart_fn,
                           '-maxh', str(opt['cube_run_time'])
                           ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)
                p.communicate()

        # Read the current_number of iterations
        current_md_steps = gmx_steps("-cp", restart_fn)

        os.chdir(cwd)

        # Upload Restart files to S3
        tar_restart_fn = prefix + "_" + Gromacs.gmx_restart_prefix + '_' + str(cycle_id) + '.tar'

        with tarfile.open(tar_restart_fn, mode='w:gz') as archive:
            archive.add(restart_fn, arcname=os.path.basename(restart_fn))
            archive.add(fn_tpr, arcname=os.path.basename(fn_tpr))

        file_restart_id = upload_file(tar_restart_fn,
                                      orion_ui_name=prefix + "_" + "Gmx_Restart" + '_' + str(cycle_id) + '.tar',
                                      tag="Gmx_Restart")

        # Upload the full gmx trajectory
        tar_trj_dir_fn = prefix + "_" + Gromacs.gmx_traj_dir_prefix + '_' + str(cycle_id) + '.tar'

        with tarfile.open(tar_trj_dir_fn, mode='w:gz') as archive:
            archive.add(outdir, arcname=os.path.basename(outdir))

        file_traj_id = upload_file(tar_trj_dir_fn,
                                   orion_ui_name=prefix + "_" + "Gmx_Trajectory" + '_' + str(cycle_id) + '.tar',
                                   tag="Gmx_Trajectory")

        record.set_value(Fields.gmx_restart, file_restart_id)
        record.set_value(Fields.trajectory, file_traj_id)

    return current_md_steps

