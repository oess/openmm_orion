#!/usr/bin/env python

from floe.api import (WorkFloe, ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes import ParallelClusterOETrajCube, ParallelMakeClusterTrajOEMols

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

clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols')

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, clusCube, clusOEMols, ofs)

ifs.success.connect(clusCube.intake)
clusCube.success.connect(clusOEMols.intake)
clusOEMols.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
