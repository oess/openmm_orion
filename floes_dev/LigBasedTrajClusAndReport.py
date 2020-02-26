#!/usr/bin/env python

from floe.api import (WorkFloe, ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes import (ConformerGatheringData,
                                       ParallelClusterOETrajCube,
                                       ParallelMakeClusterTrajOEMols,
                                       ParallelMDTrajAnalysisClusterReport,
                                       ParallelConcatenateTrajMMPBSACube,
                                       MDFloeReportCube)

job = WorkFloe("Testing Traj OEMol Clustering on a conformer")

job.description = """
Testing Ligand Clustering Floe on a conformer
#
Ex. python floes/ConfBasedTrajClustering.py --in STMD_TrajOEMol.oedb  --out STMD_LigClus.oedb
#
Parameters:
-----------
in (.oedb file): file of the MD results with Traj OEMols
#
Outputs:
--------
ofs (.oedb file): file of the MD results with Traj OEMol Clustering on a conformer.
"""

# job.uuid = "7cacc2af-cae7-4dc7-8956-fcf539861e3d"

ifs = DatasetReaderCube("ifs")

ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

confGather = ConformerGatheringData("Gathering Conformer Records")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols')
trajMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube')
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
report = MDFloeReportCube("report", title="Floe Report")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, confGather,
              clusCube, trajMMPBSA, clusOEMols, report_gen, report,
              ofs)

ifs.success.connect(confGather.intake)
confGather.success.connect(clusCube.intake)
clusCube.success.connect(trajMMPBSA.intake)
trajMMPBSA.success.connect(clusOEMols.intake)
clusOEMols.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
