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

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ParallelTrajToOEMolCube,
                                                      ParallelTrajInteractionEnergyCube,
                                                      ParallelTrajPBSACube)

job = WorkFloe("Calculate Trajectory Protein-Ligand MMPBSA Energies from Short Trajectory MD")

job.description = """
Analyse the trajectory from Protein-Ligand MD in terms of MMPBSA interaction energies between the
between the ligand and the protein after fitting the trajectory based on active site C_alphas.
"""

job.uuid = "d00de553-5f78-4496-ae96-9c8adc527f53"

# Input MD Dataset
iMDInput = DatasetReaderCube("MDInputReader", title="MD Input Reader")
iMDInput.promote_parameter("data_in", promoted_name="in",
                           title="MD Input Dataset", description="MD Input Dataset")

# This Cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube", title="Trajectory To OEMols")
IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube", title="MM Energies")
PBSACube = ParallelTrajPBSACube("TrajPBSACube", title="PBSA Energies")

trajproc_group = ParallelCubeGroup(cubes=[trajCube, IntECube, PBSACube])
job.add_group(trajproc_group)

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

job.add_cubes(iMDInput, coll_open,
              trajCube, IntECube, PBSACube,
              coll_close, check_rec,  ofs, fail)

# Success Connections
iMDInput.success.connect(coll_open.intake)
coll_open.success.connect(trajCube.intake)
trajCube.success.connect(IntECube.intake)
IntECube.success.connect(PBSACube.intake)
PBSACube.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail Connections
coll_open.failure.connect(check_rec.fail_in)
trajCube.failure.connect(check_rec.fail_in)
IntECube.failure.connect(check_rec.fail_in)
PBSACube.failure.connect(check_rec.fail_in)
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
