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

from floe.api import (WorkFloe,
                      ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.System.cubes import (ParallelRecordSizeCheck,
                                  CollectionSetting)

from MDOrion.TrjAnalysis.cubes_hintAnalysis import (ParallelComparePoseBintsToTrajBints)

job = WorkFloe("Assess Pose Stability from Binding Interactions")

job.description = """
Compare the initial pose of the ligand to the
ligand trajectory in terms of the OEHint binding interactions.
"""

#job.uuid = "d00de553-5f78-4496-ae96-9c8adc527f53"

# Input MD Dataset
AnlysInput = DatasetReaderCube("AnalysisInputReader", title="Analysis Input Reader")
AnlysInput.promote_parameter("data_in", promoted_name="in",
                           title="Analysis Input Dataset", description="Analysis Input Dataset")

# This Cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

trajBints = ParallelComparePoseBintsToTrajBints("TrajBintsCube", title="Trajectory Binding Interactions")

# This Cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success", title="Record Check Success")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(AnlysInput, coll_open,
              trajBints,
              coll_close, check_rec,  ofs, fail)

# Success Connections
AnlysInput.success.connect(coll_open.intake)
coll_open.success.connect(trajBints.intake)
trajBints.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail Connections
coll_open.failure.connect(check_rec.fail_in)
trajBints.failure.connect(check_rec.fail_in)
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
