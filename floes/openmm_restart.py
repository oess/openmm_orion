from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import OEBSinkCube

job = WorkFloe("RestartOpenMMSimulation")

job.description = """
**Restart an OpenMM Simulation from a saved state**

Ex. `python floes/openmm_restart.py --complex output/9PC1X-simulation.oeb.gz --steps 10000`

Parameters:
-----------
complex (file): OEB file of the simulated protein:ligand complex containing a saved state

Optional:
--------
steps (int): Number of MD steps to equilibrate the complex (default: 50,000)

Outputs:
--------
ofs: Outputs to a <idtag>-restart.oeb.gz file
"""

job.classification = [["Testing", "OpenMM"], ["Testing", "Simulation"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

md = OpenMMSimulation('md')
md.promote_parameter('steps', promoted_name='steps')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix="restart")

job.add_cubes(ifs, md, ofs)
ifs.success.connect(md.intake)
md.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
