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

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from floes.SubfloeFunctions import setup_MD_startup

from MDOrion.System.cubes import MDComponentCube

from MDOrion.System.cubes import ParallelSolvationCube

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelRecordSizeCheck)

job = WorkFloe('Solvate and Run MD',
               title='Solvate and Run MD')

job.description = open(path.join(path.dirname(__file__), 'PlainMD_desc.rst'), 'r').read()
# Locally the floe can be invoked by running the terminal command:
# python floes/PlainMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['General MD']]
job.uuid = "266481fc-b257-41e9-b2f9-a92bf028b701"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="solute", title='Solute Input File',
                      description="Solute input file")

sysid = IDSettingCube("System Ids")

md_comp = MDComponentCube("MD Components")
md_comp.set_parameters(multiple_flasks=True)

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = ParallelSolvationCube("Hydration", title="Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')
solvate.set_parameters(close_solvent=True)

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

job.add_cubes(ifs, sysid, md_comp, solvate, coll_open,
              coll_close, check_rec, ofs, fail)

# Connections before setup_MD_startup subfloe
ifs.success.connect(sysid.intake)
sysid.success.connect(md_comp.intake)
sysid.failure.connect(check_rec.fail_in)
md_comp.success.connect(solvate.intake)
md_comp.failure.connect(check_rec.fail_in)
solvate.success.connect(coll_open.intake)
solvate.failure.connect(check_rec.fail_in)
coll_open.failure.connect(check_rec.fail_in)

# Call subfloe function to start up the MD, equilibrate, and do the production run
MD_startup_options = {}
MD_startup_options['Prod_Default_Time_ns'] = 2.0
MD_startup_options['Prod_Default_Traj_Intvl_ns'] = 0.004
MD_outcube = setup_MD_startup(job, coll_open, check_rec, MD_startup_options)

# Connections after setup_MD_startup subfloe
MD_outcube.success.connect(coll_close.intake)
MD_outcube.failure.connect(check_rec.fail_in)
coll_close.success.connect(check_rec.intake)
coll_close.failure.connect(check_rec.fail_in)
check_rec.success.connect(ofs.intake)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
