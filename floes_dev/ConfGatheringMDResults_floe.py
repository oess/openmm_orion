#!/usr/bin/env python

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_trajProcessing import ConformerGatheringData

job = WorkFloe("Test Gathering Traj OEMol Confs")

job.description = """
Test Gathering Traj OEMol Confs
#
Ex. python floes_dev/ConfGatheringMDResults_floe.py --in  STMD_TrajOEMol.oedb
--out STMD_ConfsGathered.oedb
#
Parameters:
-----------
in (.oedb file): file of the MD results with one record per conf 
#
Outputs:
--------
ofs (.oedb file): file of the MD results with one multiconf record per ligand.
"""

# job.uuid = "7cacc2af-cae7-4dc7-8956-fcf539861e3d"

ifs = DatasetReaderCube("ifs")

ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

confGather = ConformerGatheringData("Gathering Conformer Records")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

job.add_cubes(ifs, confGather, ofs)

ifs.success.connect(confGather.intake)
confGather.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
