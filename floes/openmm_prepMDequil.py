from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from OpenMMCubes.cubesPrepMD import OpenMMequilCube
from LigPrepCubes.cubes import OEBSinkCube

job = WorkFloe("PrepMDequilibrate")

job.description = """
**Equilibrate a warmed up OpenMM-ready solvated complex**

Ex. `data='examples/data'; python floes/openmm_prepMDequil.py --complex $data/9PC1X-complex.oeb.gz --picosec 10`

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
picosec (int): Number of picoseconds to warm up the complex.
temperature (decimal): target final temperature after warming.
restraintWt (decimal): strength in kcal/mol/ang^2 for xyz atom restraints.
snapFreq (int): frequency (in picoseconds) for taking snapshots.

Outputs:
--------
ofs: Outputs to a <idtag>-equil.oeb.gz file
"""

job.classification = [['PrepMDequilibrate']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

equil = OpenMMequilCube('equil')
equil.promote_parameter('picosec', promoted_name='picosec')
equil.promote_parameter('restraintType', promoted_name='restraintType')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='equil')

job.add_cubes(ifs, equil, ofs)
ifs.success.connect(equil.intake)
equil.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
