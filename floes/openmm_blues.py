from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from OpenMMCubes.blues import BluesNCMC
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

job = WorkFloe("RunOpenMMBLUES")

job.description = """
**Using a prepared protine:ligand complex, run a short equil. MD simulation,
and then run BLUES**

This floe will do the following in each cube:
  (1) ifs: Read in the protein:ligand complex file (9PC1X-complex.oeb.gz),
  (2) blues: Run 25,000 BLUES steps.
  (3) ofs: Write out the OEMOl of the simulated complex to a <idtag>-blues.oeb.gz

Ex. `python floes/openmm_mdblues.py --complex examples/data/9PC1X-simulation.oeb.gz --mdsteps 5000 --ncsteps 25 --nciter 10`

Parameters:
-----------
complex (file): OEB file of the simulated protein:ligand complex

Optional:
--------
mdsteps (int): Number of MD steps to use during BLUES run (default: 25,000)
ncsteps (int): Number of NCMC steps to use during BLUES run (default: 25)
nciter (int): Number of iterations to perform during NCMC steps (default: 10)

Outputs:
--------
ofs: Outputs to a <idtag>-blues.oeb.gz file
"""

job.classification = [["OpenMM", "Simulation", "Enhanced Sampling"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")


blues = BluesNCMC('blues')
blues.promote_parameter('mdsteps', promoted_name='mdsteps')
blues.promote_parameter('ncsteps', promoted_name='ncsteps')
blues.promote_parameter('nciter', promoted_name='nciter')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='blues')

job.add_cubes(ifs, mblues, ofs)
ifs.success.connect(md.intake)
md.success.connect(blues.intake)
blues.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
