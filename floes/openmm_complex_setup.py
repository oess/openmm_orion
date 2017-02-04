from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

job = WorkFloe("SetupOpenMMComplex")

job.description = """
**Set up an OpenMM complex for simulation**

This floe will do the following in each cube:
  (1) ifs: Read in the ligand file (toluene.pdb),
  (2) idtag: Add an idtag from the molecule's title or use a random 6 character string.
  (3a) smirff: Parameterize the molecule with the ffxml file (smirff99Frosst.ffxml)
        Generate the ParmEd Structure and attach it to the OEMol.
  (3b) smirff_out: Write out the OEMol to a <idtag>-smirff-oeb.gz
  (4) complex_setup: Paramterize the protein (T4-protein.pdb) and merge with the molecule Structure,
        Using PDBFixer: add missing atoms, add hydrogens given a pH, and solvate with TIP3P.
        Attach tagged data containing the <idtag>, <Structure> and <System>
  (5) ofs: Write out the OEMOl of the complex to a <idtag>-complex.oeb.gz

Parameters:
-----------
ligand (ifs): expected to be a ligand already docked to the protein structure.
protein: Assumed to be taken from a 'pre-prepared' protein structure.
ffxml: The smirff99Frosst.ffxml file.

Outputs:
--------
smirff_out: Outputs a <idtag>-smirff.oeb.gz file containing the
OpenMM System and ParmEd Structure of the ligand, packaged with the OEMol.

ofs: Outputs a <idtag>-complex.oeb.gz file containing the
OpenMM System and ParmEd Structure of the protein:ligand complex,
packaged with the OEMol.
"""

job.classification = [
    ["OpenMM", "Protein-Ligand Complex Setup"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="docked ligands")

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")
smirff_out = OEBSinkCube('smirff_out')
smirff_out.set_parameters(suffix='smirff')

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_forcefield')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_forcefield')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='complex')

job.add_cubes(ifs, idtag, smirff, smirff_out, complex_setup, ofs)
ifs.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)

smirff.success.connect(smirff_out.intake)
smirff.success.connect(complex_setup.intake)

complex_setup.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
