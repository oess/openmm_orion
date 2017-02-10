from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SmilesSimulation")

job.description = """
This floe will do the following in each cube:
  (1) ifs: Read in the SMILES from file (test_smiles.ism)
  (2) omega: Generate multiconformer molecules
  (3) fred: Dock the MCMol to a prepared receptor (test-receptor.oeb.gz)
        Emit top scoring pose and attach score as SDData.
  (4) idtag: Add an idtag from the molecule's title or use a random 6 character string.
  (5) smirff: Parameterize the molecule with the ffxml file (smirff99Frosst.ffxml)
        Generate the ParmEd Structure and attach it to the OEMol.
  (6) complex_setup: Paramterize the protein (receptor-fixed.pdb) and merge with the molecule Structure,
        Using PDBFixer: add missing atoms, add hydrogens given a pH, and solvate with TIP3P.
        Attach tagged data containing the <idtag>, <Structure> and <System>.
  (7) md: Minimize the complex and run 50,000 steps of MD using the prepared complex and report every 1000 steps.
      Reporters: Progress of the simulation, state data for energies, checkpoints, DCD and h5.
      Attach tagged data containing the <idtag>, <Structure>, <System>, <State>, and <logfile>.
  (8) ofs: Write out the OEMOl of the simulated complex to a <idtag>-simulation.oeb.gz

Ex. `python floes/smiles_setup-md.py --ligand examples/data/test_smiles.ism --receptor examples/data/test-receptor.oeb.gz --protein examples/data/receptor-fixed.pdb --steps 10000`

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
Writes out files at each cube stage:
<idtag>-docked.oeb.gz
<idtag>-smirff.oeb.gz
<idtag>-complex.oeb.gz
ofs: Outputs to a <idtag>-simulation.oeb.gz
"""

job.classification = [["Simulation", "OpenMM", "Testing", "Ligand Preparation", "Complex Setup"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="File containing SMILES")

omega = OEOmegaConfGen('omega')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')
fred_out = OEBSinkCube('fred_out')
fred_out.set_parameters(suffix='docked')

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")
smirff_out = OEBSinkCube('smirff_out')
smirff_out.set_parameters(suffix='smirff')

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein', description="PDB of protein structure")
complex_setup.promote_parameter('pH', promoted_name='pH')
complex_setup.promote_parameter('solvent_padding', promoted_name='solvent_padding')
complex_setup.promote_parameter('salt_concentration', promoted_name='salt_conc')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_ff')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_ff')
complex_out = OEBSinkCube('complex_out')
complex_out.set_parameters(suffix='complex')

md = OpenMMSimulation('md')
md.promote_parameter('steps', promoted_name='steps')
ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='simulation')

cubes = [ifs, omega, fred, idtag, smirff, smirff_out,
        complex_setup, complex_out, md, ofs]

job.add_cubes(*cubes)

ifs.success.connect(omega.intake)
omega.success.connect(fred.intake)
fred.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)

smirff.success.connect(smirff_out.intake)
smirff.success.connect(complex_setup.intake)

complex_setup.success.connect(complex_out.intake)
complex_setup.success.connect(md.intake)

md.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
