from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SmilesLigPrep")

job.description = """
This floe will do the following in each cube:
  (1) ifs: Read in the SMILES from file (test_smiles.ism)
  (2) omega: Generate multiconformer molecules
  (3) fred: Dock the MCMol to a prepared receptor (test-receptor.oeb.gz)
        Emit top scoring pose and attach score as SDData.
  (4) idtag: Add an idtag from the molecule's title or use a random 6 character string.
  (5) smirff: Parameterize the molecule with the ffxml file (smirff99Frosst.ffxml)
        Generate the ParmEd Structure and attach it to the OEMol.
  (6) ofs: Write out the OEMOl of the complex to a <idtag>-complex.oeb.gz

Ex. `python floes/smiles_ligprep.py --ligand examples/data/test_smiles.ism --receptor examples/data/test-receptor.oeb.gz`

Parameters:
-----------
ligand (file): .ISM file containing SMILE strings
receptor (file): OEB of a receptor prepared for docking.

*Optionals:
-----------
molecule_forcefield (file): Smarty parsable FFXML file containining parameters for the molecule (default: smirff99Frosst.ffxml)

Outputs:
--------
ofs: Outputs a <idtag>-smirff.oeb.gz file containing: <idtag>, <Structure> and <System>.
attached to the OEMol of the ligand as generic data.
"""

job.classification = [["Ligand Preparation"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="File containing SMILES")

omega = OEOmegaConfGen('omega')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='smirff')

job.add_cubes(ifs, omega, fred, idtag, smirff, ofs)
ifs.success.connect(omega.intake)
omega.success.connect(fred.intake)
fred.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)
smirff.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
