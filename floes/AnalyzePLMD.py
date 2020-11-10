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

from os import path

from floe.api import (WorkFloe,
                      ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.System.cubes import (ParallelRecordSizeCheck)

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ParallelTrajToOEMolCube,
                                                      ParallelTrajInteractionEnergyCube,
                                                      ParallelTrajPBSACube,
                                                      ConformerGatheringData,
                                                      ParallelConfTrajsToLigTraj,
                                                      ParallelConcatenateTrajMMPBSACube)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelMakeClusterTrajOEMols,
                                                       ParallelMDTrajAnalysisClusterReport,
                                                       ParallelClusterPopAnalysis,
                                                       ParallelTrajAnalysisReportDataset,
                                                       MDFloeReportCube)

from MDOrion.System.cubes import CollectionSetting

job = WorkFloe('Analyze Protein-Ligand MD',
               title='Analyze Protein-Ligand MD')

job.description = open(path.join(path.dirname(__file__), 'AnalyzePLMD_desc.rst'), 'r').read()

job.classification = [['Specialized MD']]
job.uuid = "7438db4d-30b1-478c-afc0-e921f0336c78"
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
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

confGather = ConformerGatheringData("Gathering Conformer Records",  title="Gathering Conformer Records")
catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj", title="Combine Pose Trajectories")
catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube', title="Concatenate MMPBSA Energies")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube", title="Clustering")
clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis',  title="Clustering Analysis")
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols', title="Per-Cluster Analysis")
prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset', title="Analysis Report")
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport", title="Relevant Output Extraction")

analysis_group = ParallelCubeGroup(cubes=[catLigTraj, catLigMMPBSA, clusCube, clusPop,
                                          clusOEMols, prepDataset, report_gen])
job.add_group(analysis_group)

report = MDFloeReportCube("report", title="Floe Report")

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
              trajCube, IntECube, PBSACube, confGather,
              catLigTraj, catLigMMPBSA, clusCube, clusPop, clusOEMols,
              prepDataset, report_gen, report,
              coll_close, check_rec,  ofs, fail)

# Success Connections
iMDInput.success.connect(coll_open.intake)
coll_open.success.connect(trajCube.intake)
trajCube.success.connect(IntECube.intake)
IntECube.success.connect(PBSACube.intake)
PBSACube.success.connect(confGather.intake)
confGather.success.connect(catLigTraj.intake)
catLigTraj.success.connect(catLigMMPBSA.intake)
catLigMMPBSA.success.connect(clusCube.intake)
clusCube.success.connect(clusPop.intake)
clusPop.success.connect(clusOEMols.intake)
clusOEMols.success.connect(prepDataset.intake)
prepDataset.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail Connections
coll_open.failure.connect(check_rec.fail_in)
trajCube.failure.connect(check_rec.fail_in)
IntECube.failure.connect(check_rec.fail_in)
PBSACube.failure.connect(check_rec.fail_in)
confGather.failure.connect(check_rec.fail_in)
catLigTraj.failure.connect(check_rec.fail_in)
catLigMMPBSA.failure.connect(check_rec.fail_in)
clusCube.failure.connect(check_rec.fail_in)
clusPop.failure.connect(check_rec.fail_in)
clusOEMols.failure.connect(check_rec.fail_in)
prepDataset.failure.connect(check_rec.fail_in)
report_gen.failure.connect(check_rec.fail_in)
report.failure.connect(check_rec.fail_in)
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
