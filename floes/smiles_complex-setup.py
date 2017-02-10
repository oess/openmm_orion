from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SmilesComplexPrep")

job.description = """
This floe will do the following in each cube:
  (1) ifs: Read in the SMILES from file (test_smiles.ism)
  (2) omega: Generate multiconformer molecules
  (3) fred: Dock the MCMol to a prepared receptor (test-receptor.oeb.gz)
        Emit top scoring pose and attach score as SDData.
  (4) idtag: Add an idtag from the molecule's title or use a random 6 character string.
  (5) smirff: Parameterize the molecule with the ffxml file (smirff99Frosst.ffxml)
        Generate the ParmEd Structure and attach it to the OEMol.
  (6) complex_setup: Paramterize the protein (T4-protein.pdb) and merge with the molecule Structure,
        Using PDBFixer: add missing atoms, add hydrogens given a pH, and solvate with TIP3P.
        Attach tagged data containing the <idtag>, <Structure> and <System>.
  (7) ofs: Write out the OEMOl of the complex to a <idtag>-complex.oeb.gz

Ex. `python floes/smiles_complex-setup.py --ligand examples/data/test_smiles.ism --receptor examples/data/test-receptor.oeb.gz --protein examples/data/receptor-fixed.pdb`

Parameters:
-----------
ligand (file): .ISM file containing SMILE strings
receptor (file): OEB of a receptor prepared for docking.
protein (file): PDB of the prepared protein structure.

*Optionals:
-----------
pH (float): Solvent pH used to select protein protonation states (default: 7.0)
solvent_padding (float): Padding around protein for solvent box (default: 10 angstroms)
salt_concentration (float): Salt concentration (default: 50 millimolar)
molecule_forcefield (file): Smarty parsable FFXML file containining parameters for the molecule (default: smirff99Frosst.ffxml)
protein_forcefield (file): XML file containing forcefield parameters for protein (default: amber99sbildn.xml)
solvent_forcefield (file): XML file containing forcefield parameter for solvent (default: tip3p.xml)

Outputs:
--------
ofs: Outputs a <idtag>-complex.oeb.gz file containing: <idtag>, <Structure> and <System>.
attached to the OEMol of the protein:ligand complex as generic data.
"""

job.classification = [
    ["Testing", "Complex Setup"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="File containing SMILES")

omega = OEOmegaConfGen('omega')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein', description="PDB of protein structure")
complex_setup.promote_parameter('pH', promoted_name='pH')
complex_setup.promote_parameter('solvent_padding', promoted_name='solvent_padding')
complex_setup.promote_parameter('salt_concentration', promoted_name='salt_conc')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_ff')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_ff')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='complex')

job.add_cubes(ifs, omega, fred, idtag, smirff, complex_setup, ofs)
ifs.success.connect(omega.intake)
omega.success.connect(fred.intake)
fred.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)
smirff.success.connect(complex_setup.intake)
complex_setup.success.connect(ofs.intake)
if __name__ == "__main__":
    job.run()
