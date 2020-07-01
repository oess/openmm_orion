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

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelMDTrajAnalysisClusterReport,
                                                       MDFloeReportCube)

job = WorkFloe('Short Trajectory MD cluster and report',
               title='Short Trajectory MD cluster and report')

job.description = """
An MD production run which has been analyzed in terms of interactions between the
ligand and the active site is clustered in terms of ligand RMSD after fitting the trajectory
based on active site C_alphas.

Required Input Parameters:
--------------------------
Dataset result of an MD Short Trajectory run with traj OEMols and MMPBSA interactions

Outputs:
--------
floe report: html page of the Analysis of each ligand.
out (.oedb file): file of the Analysis results for all ligands.
"""

job.classification = [['Analysis']]
job.uuid = "4de8826b-8902-4037-9f8e-5c1a7d8a793c"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

report = MDFloeReportCube("report", title="Floe Report")

job.add_cubes(ifs, clusCube, report_gen, report, ofs, fail)

ifs.success.connect(clusCube.intake)
clusCube.success.connect(report_gen.intake)
clusCube.failure.connect(fail.intake)
report_gen.success.connect(report.intake)
report_gen.failure.connect(fail.intake)
report.failure.connect(fail.intake)
report.success.connect(ofs.intake)
report.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
