#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelClusterPopAnalysis,
                                                       ParallelMakeClusterTrajOEMols)

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ConformerGatheringData,
                                                      ParallelConfTrajsToLigTraj,
                                                      ParallelConcatenateTrajMMPBSACube)

from MDOrion.TrjAnalysis.cubes_hintAnalysis import (ParallelBintScoreInitialPoseAndTrajectory)

job = WorkFloe("Test from combining confs Traj OEMols through Clustering")

job.description = """
Testing from aggregating conf trajs (into one ligand traj OEMol) through clustering.
The input dataset is an .oedb file of the unaggregated confs MD results with Traj OEMols + IntE + PBSA
"""

# job.uuid = "7cacc2af-cae7-4dc7-8956-fcf539861e3d"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

confGather = ConformerGatheringData("Gathering Conformer Records", title="Gathering Conformer Records")
catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj", title="Combine Pose Trajectories")
catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube', title="Concatenate MMPBSA Energies")
trajBints = ParallelBintScoreInitialPoseAndTrajectory("TrajBintsCube", title="Trajectory Binding Interactions")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube", title="Clustering")
clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis', title="Clustering Analysis")
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols', title="Per-Cluster Analysis")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs,
              confGather, catLigTraj, catLigMMPBSA, trajBints, clusCube, clusPop, clusOEMols,
              ofs, fail)

ifs.success.connect(confGather.intake)
confGather.success.connect(catLigTraj.intake)
catLigTraj.success.connect(catLigMMPBSA.intake)
catLigMMPBSA.success.connect(trajBints.intake)
trajBints.success.connect(clusCube.intake)
clusCube.success.connect(clusPop.intake)
clusPop.success.connect(clusOEMols.intake)
clusOEMols.success.connect(ofs.intake)
clusOEMols.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
