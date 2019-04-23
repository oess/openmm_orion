#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
from MDOrion.TrjAnalysis.cubes import TrajInteractionEnergyCube

job = WorkFloe("Testing Traj Protein-Ligand Interaction Energies")

job.description = """
Testing Protein-Ligand Interaction Energies Floe
#
Ex. python floes/up.py --in  STMD_TrajOEMol.oedb
--out STMD_TrajIntE.oedb
#
Parameters:
-----------
in (.oedb file): file of the MD results with Traj OEMols
#
Outputs:
--------
ofs (.oedb file): file of the MD results with Interaction Energy results.
"""

ifs = DatasetReaderCube("ifs")

ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

scube = TrajInteractionEnergyCube("TrajInteractionEnergyCube")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, scube, ofs)

ifs.success.connect(scube.intake)
scube.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
