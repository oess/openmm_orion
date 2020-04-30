#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes import ParallelTrajToOEMolCube

from MDOrion.TrjAnalysis.cubes import ParallelTrajInteractionEnergyCube

from MDOrion.TrjAnalysis.cubes import ParallelTrajPBSACube

job = WorkFloe("Trajectory Interaction Energies from Short Trajectory MD")

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

job.uuid = "d00de553-5f78-4496-ae96-9c8adc527f53"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube")
trajIntE = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube")
trajPBSA = ParallelTrajPBSACube("TrajPBSACube")

job.add_cubes(ifs, trajCube, trajIntE, trajPBSA, ofs)

ifs.success.connect(trajCube.intake)
trajCube.success.connect(trajIntE.intake)
trajIntE.success.connect(trajPBSA.intake)
trajPBSA.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
