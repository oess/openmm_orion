#!/usr/bin/env python

from floe.api import (WorkFloe, ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (
                                    ParallelClusterPopAnalysis,
                                    ParallelMDTrajAnalysisClusterReport,
                                    MDFloeReportCube)

job = WorkFloe("Starting with STMD results through clustering, just generate the Analysis Floe report")

job.description = """
Starting with STMD results through clustering, just generate the Analysis Floe report.
The input dataset is an .oedb file of the STMD results through clustering
"""

# job.uuid = "7cacc2af-cae7-4dc7-8956-fcf539861e3d"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

popanlys = ParallelClusterPopAnalysis('ClusterPopAnalysis')
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
report = MDFloeReportCube("report", title="Floe Report")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs,
              popanlys,
              report_gen, report,
              ofs, fail)

ifs.success.connect(popanlys.intake)
popanlys.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report_gen.failure.connect(fail.intake)
report.success.connect(ofs.intake)
report.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
