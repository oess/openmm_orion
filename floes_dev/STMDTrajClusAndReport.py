#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube


from MDOrion.TrjAnalysis.cubes_trajProcessing import (ParallelTrajInteractionEnergyCube,
                                                      ParallelTrajPBSACube)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelMDTrajAnalysisClusterReport,
                                                       MDFloeReportCube)

job = WorkFloe("Analysing Trajectory from Short Trajectory MD")

job.description = """
Analyse the trajectory from Short Trajectory MD in terms of interactions between the
ligand and the active site and in terms of ligand RMSD after fitting the trajectory
based on active site C_alphas.

Required Input Parameters:
--------------------------
in: Collection of OERecords (one per ligand) of Short Trajectory MD results.

Outputs:
--------
floe report: html page of the Analysis for each ligand.
out (.oedb file): file of the Analysis results for all ligands.
"""

job.uuid = "673da3f5-3612-4a0d-a964-688e93b1e319"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

# trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube")
trajIntE = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube")
trajPBSA = ParallelTrajPBSACube("TrajPBSACube")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
molHtml = ParallelMDTrajAnalysisClusterReport("MolHtmlCube")
floeReport = MDFloeReportCube("FloeReportCube")

# job.add_cubes(ifs, trajCube, trajIntE, trajPBSA, clusCube, molHtml, ofs)
job.add_cubes(ifs, trajIntE, trajPBSA, clusCube, molHtml, floeReport, ofs)

# ifs.success.connect(trajCube.intake)
# trajCube.success.connect(trajIntE.intake)
ifs.success.connect(trajIntE.intake)
trajIntE.success.connect(trajPBSA.intake)
trajPBSA.success.connect(clusCube.intake)
clusCube.success.connect(molHtml.intake)
molHtml.success.connect(floeReport.intake)
floeReport.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
