from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubesPrepMD import OpenMMminimizeCube

job = WorkFloe("PrepMDminimize")

job.description = """
**Minimize an OpenMM-ready solvated complex**

Ex. `data='examples/data'; python floes/openmm_prepMDminimize.py --complex $data/9PC1X-complex.oeb.gz --steps 1000`

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
steps (int): Number of MD steps to equilibrate the complex (default: 50,000)

Outputs:
--------
ofs: Outputs to a <idtag>-simulation.oeb.gz file
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

minComplex = OpenMMminimizeCube('minComplex')
minComplex.promote_parameter('steps', promoted_name='steps')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, minComplex, ofs, fail)
ifs.success.connect(minComplex.intake)
minComplex.success.connect(ofs.intake)
minComplex.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
