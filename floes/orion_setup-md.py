from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("OrionSimulation")

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

Ex. `python floes/orion_setup-md.py --ligand input/test_smiles.ism --receptor input/test-receptor.oeb.gz --ffxml input/smirff99Frosst.ffxml --protein input/receptor-fixed.pdb --steps 10000`

Parameters:
-----------
ligand(ifs): .ISM file containing SMILE strings
receptor: .OEB of a receptor prepared for docking.
ffxml: The smirff99Frosst.ffxml file.
protein: prepared PDB file of the receptor

Optional:
--------
steps: Number of MD steps to equilibrate the complex (default: 50,000)


Outputs:
--------
Writes out files at each cube stage:
<idtag>-docked.oeb.gz
<idtag>-smirff.oeb.gz
<idtag>-complex.oeb.gz
ofs: Outputs to a <idtag>-simulation.oeb.gz
"""

job.classification =[["Testing", "Simulation", "OpenMM"]]
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
smirff.set_parameters(molecule_forcefield='smirff99Frosst.ffxml')
smirff_out = OEBSinkCube('smirff_out')
smirff_out.set_parameters(suffix='smirff')

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein')
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
