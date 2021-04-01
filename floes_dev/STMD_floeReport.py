#!/usr/bin/env python
from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (
                                       ParallelTrajAnalysisReportDataset,
                                       ParallelMDTrajAnalysisClusterReport,
                                       MDFloeReportCube)

job = WorkFloe("Floe Report from Analyzed Short Trajectory MD")

job.description = """
Generate a Floe Report from Analyzed and Clustered Short Trajectory MD results

Required Input Parameters:
--------------------------
in: Collection of OERecords (one per ligand) of analyzed, clustered Short Trajectory MD results.

Outputs:
--------
floe report: html page of the Analysis for each ligand.
out (.oedb file): file of the Analysis results for all ligands.
"""

job.uuid = "a21708dd-085b-4494-84db-a45d12c4dded"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset', title="Analysis Report")
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport", title="Relevant Output Extraction")
report = MDFloeReportCube("report", title="Floe Report")

job.add_cubes(ifs, prepDataset, report_gen, report, ofs)

ifs.success.connect(prepDataset.intake)
prepDataset.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
