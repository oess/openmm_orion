#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
from MDOrion.TrjAnalysis.cubes import TrajToOEMolCube
from MDOrion.TrjAnalysis.cubes import TrajInteractionEnergyCube
from MDOrion.TrjAnalysis.cubes import TrajPBSACube
#
job = WorkFloe("Trajectory Interaction Energies from Short Trajectory MD")
#
job.description = """
Analyse the trajectory from Short Trajectory MD in terms of interactions between the
ligand and the active site after fitting the trajectory based on active site C_alphas.

Required Input Parameters:
--------------------------
in: Collection of OERecords (one per ligand) of Short Trajectory MD results.

Outputs:
--------
out (.oedb file): file of the Analysis results for all ligands.
"""

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

trajCube = TrajToOEMolCube("TrajToOEMolCube")
trajIntE = TrajInteractionEnergyCube("TrajInteractionEnergyCube")
trajPBSA = TrajPBSACube("TrajPBSACube")

job.add_cubes(ifs, trajCube, trajIntE, trajPBSA, ofs)

ifs.success.connect(trajCube.intake)
trajCube.success.connect(trajIntE.intake)
trajIntE.success.connect(trajPBSA.intake)
trajPBSA.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
