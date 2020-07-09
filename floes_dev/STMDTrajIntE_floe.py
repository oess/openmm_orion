#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ParallelTrajToOEMolCube,
                                                      ParallelTrajInteractionEnergyCube,
                                                      ParallelTrajPBSACube)

job = WorkFloe("Calculate Trajectory Protein-Ligand MMPBSA Energies from Short Trajectory MD")

job.description = """
Analyse the trajectory from Short Trajectory MD in terms of MMPBSA interaction energies between the
between the ligand and the protein after fitting the trajectory based on active site C_alphas.
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
