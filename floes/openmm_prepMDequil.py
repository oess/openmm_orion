from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMComplexSetup
from OpenMMCubes.cubesPrepMD import OpenMMequilCube

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
equil.promote_parameter('picosec', promoted_name='picosec', default=10.0,
      description='Length of MD run in picoseconds')
equil.promote_parameter('snapFreq', promoted_name='snapFreq', default=0,
      description='frequency (in picoseconds) for taking snapshots')
equil.promote_parameter('restraintType', promoted_name='restraintType',
      default='NonSolventNonH', description='Which kinds of atoms get xyz restraints')
equil.promote_parameter('restraintWt', promoted_name='restraintWt',
      default=2.0, description='Restraint weight in kcal/mol/ang^2 for xyz atom restraints')
equil.promote_parameter('temperature', promoted_name='temperature', default=300.0,
      description='Temperature (Kelvin)')
equil.promote_parameter('label', promoted_name='label',
      default='_equil', description='label to add to filenames from this cube')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')

job.add_cubes(ifs, equil, ofs, fail)
ifs.success.connect(equil.intake)
equil.success.connect(ofs.intake)
equil.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
