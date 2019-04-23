#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
from MDOrion.TrjAnalysis.cubes import TrajToOEMolCube
from MDOrion.TrjAnalysis.cubes import TrajInteractionEnergyCube
from MDOrion.TrjAnalysis.cubes import TrajPBSACube
from MDOrion.TrjAnalysis.cubes import ClusterOETrajCube
from MDOrion.TrjAnalysis.cubes import MDTrajAnalysisClusterReport
#
job = WorkFloe("Analysing Trajectory from Short Trajectory MD")
#
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

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

trajCube = TrajToOEMolCube("TrajToOEMolCube")
trajIntE = TrajInteractionEnergyCube("TrajInteractionEnergyCube")
trajPBSA = TrajPBSACube("TrajPBSACube")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
reportCube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

job.add_cubes(ifs, trajCube, trajIntE, trajPBSA, clusCube, reportCube, ofs)

ifs.success.connect(trajCube.intake)
trajCube.success.connect(trajIntE.intake)
trajIntE.success.connect(trajPBSA.intake)
trajPBSA.success.connect(clusCube.intake)
clusCube.success.connect(reportCube.intake)
reportCube.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
