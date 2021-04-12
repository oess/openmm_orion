#!/usr/bin/env python

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

from os import path

from floe.api import (WorkFloe)

from floes.SubfloeFunctions import setup_MD_startup, setup_PLComplex_for_MD

from orionplatform.cubes import DatasetWriterCube

from MDOrion.System.cubes import (CollectionSetting,
                                  ParallelRecordSizeCheck)

job = WorkFloe('Solvate and Run Protein-Ligand MD', title='Solvate and Run Protein-Ligand MD')
job.description = open(path.join(path.dirname(__file__), 'ProteinLigandMD_desc.rst'), 'r').read()
job.classification = [['Specialized MD']]
job.uuid = "ae561d76-a2b6-4d89-b621-b979f1930b40"
job.tags = [tag for lists in job.classification for tag in lists]


# This Cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)
coll_open.set_parameters(write_new_collection='MD_OPLMD')

# This Cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(coll_open, coll_close, check_rec, ofs, fail)

# Call subfloe function to set up the solvated protein-ligand complex
PLComplex_for_MD_options = {}
PLComplex_for_MD_options['charge_ligands'] = True
PLComplex_outcube = setup_PLComplex_for_MD(job, check_rec, PLComplex_for_MD_options)

# Connections
PLComplex_outcube.success.connect(coll_open.intake)
PLComplex_outcube.failure.connect(check_rec.fail_in)
coll_open.failure.connect(check_rec.fail_in)

# Call subfloe function to start up the MD, equilibrate, and do the production run
MD_startup_options = {}
MD_startup_options['Prod_Default_Time_ns'] = 2.0
MD_startup_options['Prod_Default_Traj_Intvl_ns'] = 0.004
MD_outcube = setup_MD_startup(job, coll_open, check_rec, MD_startup_options)

# Connections
MD_outcube.success.connect(coll_close.intake)
MD_outcube.failure.connect(check_rec.fail_in)
coll_close.success.connect(check_rec.intake)
coll_close.failure.connect(check_rec.fail_in)
check_rec.success.connect(ofs.intake)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
