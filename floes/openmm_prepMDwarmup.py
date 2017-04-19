from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubesPrepMD import OpenMMwarmupNVTCube

job = WorkFloe("PrepMDwarmup")

job.description = """
**Warm up an OpenMM-ready solvated complex**

Ex. `data='examples/data'; python floes/openmm_prepMDwarmup.py --complex $data/9PC1X-complex.oeb.gz --picosec 10`

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
picosec (int): Number of picoseconds to warm up the complex.
temperature (decimal): target final temperature after warming.
restraintWt (decimal): strength in kcal/mol/ang^2 for xyz atom restraints.

Outputs:
--------
ofs: Outputs to a <idtag>-warmup.oeb.gz file
"""

job.classification = [['PrepMDwarmup']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

warmup = OpenMMwarmupNVTCube('warmup')
warmup.promote_parameter('picosec', promoted_name='picosec')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, warmup, ofs, fail)
ifs.success.connect(warmup.intake)
warmup.success.connect(ofs.intake)
warmup.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
