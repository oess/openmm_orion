#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
from MDOrion.TrjAnalysis.cubes import ClusterOETrajCube

job = WorkFloe("Testing Traj OEMol Clustering")

job.description = """
Testing Ligand Clustering Floe
#
Ex. python floes/up.py --in  STMD_TrajOEMol.oedb
--out STMD_LigClus.oedb
#
Parameters:
-----------
in (.oedb file): file of the MD results with Traj OEMols
#
Outputs:
--------
ofs (.oedb file): file of the MD results with Ligand Clustering results.
"""

ifs = DatasetReaderCube("ifs")

ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

scube = ClusterOETrajCube("ClusterOETrajCube")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, scube, ofs)

ifs.success.connect(scube.intake)
scube.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
