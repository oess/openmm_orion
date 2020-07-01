#!/usr/bin/env python

from floe.api import (WorkFloe, ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ConformerGatheringData,
                                                      ParallelConfTrajsToLigTraj,
                                                      ParallelConcatenateTrajMMPBSACube)

job = WorkFloe("Testing combining confs Traj OEMols to lig Traj OEMol")

job.description = """
Testing Aggregating conf trajs into one ligand traj OEMol.
The input dataset is an .oedb file of the aggregated confs MD results with Traj OEMols + IntE + PBSA
"""

# job.uuid = "7cacc2af-cae7-4dc7-8956-fcf539861e3d"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

confGather = ConformerGatheringData("Gathering Conformer Records")
ligTrajCube = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj")
ligMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube')

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs,
              confGather, ligTrajCube, ligMMPBSA,
              ofs, fail)

ifs.success.connect(confGather.intake)
confGather.success.connect(ligTrajCube.intake)
ligTrajCube.success.connect(ligMMPBSA.intake)
ligMMPBSA.success.connect(ofs.intake)
ligMMPBSA.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
