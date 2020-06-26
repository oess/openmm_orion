#!/usr/bin/env python

from floe.api import (WorkFloe)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ConformerGatheringData,
                                                      ParallelConfTrajsToLigTraj,
                                                      ParallelConcatenateTrajMMPBSACube)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelMakeClusterTrajOEMols,
                                                       ParallelMDTrajAnalysisClusterReport,
                                                       ParallelClusterPopAnalysis,
                                                       ParallelTrajAnalysisReportDataset,
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
catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj")
catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube')
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis')
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols')
prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset')
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
report = MDFloeReportCube("report", title="Floe Report")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, confGather,
              catLigTraj, catLigMMPBSA, clusCube, clusPop, clusOEMols,
              prepDataset, report_gen, report,
              ofs)

ifs.success.connect(confGather.intake)
confGather.success.connect(catLigTraj.intake)
catLigTraj.success.connect(catLigMMPBSA.intake)
catLigMMPBSA.success.connect(clusCube.intake)
clusCube.success.connect(clusPop.intake)
clusPop.success.connect(clusOEMols.intake)
clusOEMols.success.connect(prepDataset.intake)
prepDataset.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
