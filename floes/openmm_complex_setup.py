from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization, FREDDocking

job = WorkFloe("SetupOpenMMComplex")

job.description = """
**Set up an OpenMM complex for simulation.**

Read in PDB of the docked ligand, create multi-conformer molecule with OMEGA,
assign charges with QUACPAC, and parameterize the molecule with the chosen forcefields. Using PDBFixer, add missing atoms and assign protonation states, and solvate the protein:ligand complex with TIP3P.
"""

job.classification = [["Complex Setup"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="PDB of docked ligand")

charge = ChargeMCMol('charge')

lig_param = LigandParameterization('lig_param')
lig_param.promote_parameter('molecule_forcefield', promoted_name='molecule_forcefield', description='Forcefield for molecule')

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein', description="PDB of protein structure")
complex_setup.promote_parameter('pH', promoted_name='pH')
complex_setup.promote_parameter('solvent_padding', promoted_name='solvent_padding')
complex_setup.promote_parameter('salt_concentration', promoted_name='salt_conc')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_ff')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_ff')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')

job.add_cubes(ifs, charge, lig_param, complex_setup, ofs, fail)
ifs.success.connect(charge.intake)
charge.success.connect(lig_param.intake)
lig_param.success.connect(complex_setup.intake)
complex_setup.success.connect(ofs.intake)
complex_setup.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
