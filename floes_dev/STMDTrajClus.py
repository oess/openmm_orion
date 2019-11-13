#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes import ParallelClusterOETrajCube

from MDOrion.TrjAnalysis.cubes import ParallelMDTrajAnalysisClusterReport

from MDOrion.TrjAnalysis.cubes import MDFloeReportCube

job = WorkFloe("Analysing Trajectory from Short Trajectory MD")

job.description = """
Cluster the trajectory from Short Trajectory MD in terms of ligand RMSD after fitting the trajectory
based on active site C_alphas.

Required Input Parameters:
--------------------------
in: Collection of OERecords (one per ligand) of Short Trajectory MD results.

Outputs:
--------
out (.oedb file): file of the Analysis results for all ligands.
"""

job.uuid = "43f33e3f-0240-4e34-9b8b-da4d5796052a"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
molHtml = ParallelMDTrajAnalysisClusterReport("MolHtmlCube")
floeReport = MDFloeReportCube("FloeReportCube")

job.add_cubes(ifs, clusCube, molHtml, floeReport, ofs)

ifs.success.connect(clusCube.intake)
clusCube.success.connect(molHtml.intake)
molHtml.success.connect(floeReport.intake)
floeReport.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
