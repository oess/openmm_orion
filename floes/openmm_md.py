from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization, FREDDocking

job = WorkFloe("RunOpenMMSimulation")

job.description = """
**Run an OpenMM Simulation (default 500K steps / 1ns)**

Simulation cube will minimize the protein:ligand complex or restart from a Saved state.
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

md = OpenMMSimulation('md')
md.promote_parameter('steps', promoted_name='steps')


ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')

job.add_cubes(ifs, md, ofs, fail)
ifs.success.connect(md.intake)
md.success.connect(ofs.intake)
md.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
