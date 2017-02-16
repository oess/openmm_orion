from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubesPrepMD import OpenMMmakePLmaskCube

job = WorkFloe("Make Protein-Ligand Mask")

job.description = """
**Generate the Protein-Ligand Mask used for MD restraints**

Ex. `data='examples/data'; python floes/openmm_makePLmask.py --complex $data/9PC1X-complex.oeb.gz --activSiteResNums ''`

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
activeSiteResNums (file): residue numbers of active site protein residues

Outputs:
--------
ofs: Outputs to a <idtag>-PLmask.oeb.gz file
"""

job.classification = [['makePLmask']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the solvated protein:ligand complex")

PLmask = OpenMMmakePLmaskCube('PLmask')
PLmask.promote_parameter("ActSiteResNumSDTag", promoted_name="ActSiteResNumSDTag", description='whitespace delimited list of integers corresponding to residue numbers')

ofs = OEMolOStreamCube('ofs')
ofs.set_parameters(data_out="complex_PLmask.oeb.gz")

job.add_cubes(ifs, PLmask, ofs)
ifs.success.connect(PLmask.intake)
PLmask.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
