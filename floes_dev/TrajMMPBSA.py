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

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes import (ParallelTrajToOEMolCube,
                                       ParallelTrajInteractionEnergyCube,
                                       ParallelTrajPBSACube,
                                       ParallelClusterOETrajCube,
                                       ParallelMDTrajAnalysisClusterReport,
                                       MDFloeReportCube)

job = WorkFloe('Trajectory MMPBSA from Traj OEMOls',
               title='Trajectory MMPBSA from Traj OEMOls')

job.description = """
An MMPBSA analysis is carried out on trajectory OEMols for protein, ligand
and possibly waters.

Required Input Parameters:
--------------------------
.oedb of records containing protein, ligand and possibly water
trajectory OEMols from an MD Short Trajectory run.

Outputs:
--------
out (.oedb file): file of the Analysis results for all ligands.
"""

job.classification = [['Analysis']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube")
PBSACube = ParallelTrajPBSACube("TrajPBSACube")

report = MDFloeReportCube("report", title="Floe Report")

job.add_cubes(ifs, IntECube, PBSACube, ofs, fail)

ifs.success.connect(IntECube.intake)
IntECube.success.connect(PBSACube.intake)
IntECube.failure.connect(fail.intake)
PBSACube.success.connect(ofs.intake)
PBSACube.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
