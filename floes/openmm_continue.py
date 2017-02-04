from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import OEBSinkCube

job = WorkFloe("RestartOpenMMSimulation")

job.description = """
**Restart an OpenMM Simulation**

Check out the awesome stuff at the [OpenMM website](http://openmm.org)
"""

job.classification = [
    ["OpenEye", "OpenMM"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ifs", description="complex file")

md_sim = OpenMMSimulation('md_sim')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix="restart")

job.add_cubes(ifs, md_sim, ofs)
ifs.success.connect(md_sim.intake)
md_sim.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
