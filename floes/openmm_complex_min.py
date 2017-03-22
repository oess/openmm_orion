from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization, FREDDocking

job = WorkFloe("MinimizeComplex")

job.description = """
**Generate a TIP3P solvated protein:ligand system, minimize, and run a short 1ps MD trajectory (500 steps).**

*Note: Floe expects parameterized molecules (i.e. OEMols attached with a parmed.Structure).*

PDBFixer is used to add missing atoms, assign protonation states, and solvate the protein:ligand complex with TIP3P. By default the protein is parameterized using amber99sbildn forcefield parameters.
"""

job.classification = [["Complex Setup", "Minimization"],]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="docked ligands")

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_forcefield')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_forcefield')

md = OpenMMSimulation('md')
md.set_parameters(steps=500) #1ps

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')

job.add_cubes(ifs, complex_setup, md, ofs, fail)
ifs.success.connect(complex_setup.intake)
complex_setup.success.connect(md.intake)
md.success.connect(ofs.intake)
md.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
