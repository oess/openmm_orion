#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_trajProcessing import ParallelTrajPBSACube

job = WorkFloe("Testing Traj PBSA Energies")

job.description = """
Testing PBSALigand Interaction Energies Floe
#
Ex. python floes/up.py --in  STMD_TrajIntE.oedb
--out STMD_TrajPBSA.oedb
#
Parameters:
-----------
in (.oedb file): file of the Interaction Energy results with Traj OEMols
#
Outputs:
--------
ofs (.oedb file): file of the MD results with PBSA results.
"""

job.uuid = "c25a9cdb-9f84-4cfe-828e-64300a4f8c64"

ifs = DatasetReaderCube("ifs")

ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

scube = ParallelTrajPBSACube("TrajPBSACube")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, scube, ofs)

ifs.success.connect(scube.intake)
scube.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
