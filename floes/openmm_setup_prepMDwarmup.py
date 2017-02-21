from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube
from OpenMMCubes.cubesPrepMD import OpenMMmakePLmaskCube, OpenMMminimizeCube, OpenMMwarmupNVTCube, OpenMMequilCube

job = WorkFloe("Setup,Min,Warmup")

job.description = """
**Set up an OpenMM complex then run MD a la Merck Frosst**

This floe will do the following in each cube:
  (1) ifs: Read in the ligand file (toluene.pdb),
  (2) idtag: Add an idtag from the molecule's title or use a random 6 character string.
  (3a) smirff: Parameterize the molecule with the ffxml file (smirff99Frosst.ffxml)
        Generate the ParmEd Structure and attach it to the OEMol.
  (4) complex_setup: Parameterize the protein (T4-protein.pdb) and merge with the molecule Structure,
        Using PDBFixer: add missing atoms, add hydrogens given a pH, and solvate with TIP3P.
        Attach tagged data containing the <idtag>, <Structure> and <System>
  (5) Run restrained minimization and warmup similar to the Merck Frosst protocol.
  (6) ofs: Write out the OEMOl of the complex to a <idtag>-complex.oeb.gz

Ex. `python floes/openmm_setup_warmup.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb`

Parameters:
-----------
ligand (file): PDB file of ligand posed in the protein active site.
protein (file): PDB file of the protein structure, *assumed to be `pre-prepared`*

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
ofs: Outputs a <idtag>-warmup.oeb.gz file containing the
OpenMM System, State, and ParmEd Structure of the protein:ligand complex,
packaged with the OEMol.
"""

job.classification = [['Complex Setup', 'PrepMDwarmup']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="PDB of docked ligand")

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein', description="PDB of protein structure")
complex_setup.promote_parameter('pH', promoted_name='pH')
complex_setup.promote_parameter('solvent_padding', promoted_name='solvent_padding')
complex_setup.promote_parameter('salt_concentration', promoted_name='salt_conc', default=100)
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_ff')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_ff')

PLmask = OpenMMmakePLmaskCube('PLmask')
PLmask.promote_parameter("ActSiteResNumSDTag", promoted_name="ActSiteResNumSDTag", description='whitespace delimited list of integers corresponding to residue numbers')

minComplex = OpenMMminimizeCube('minComplex')
minComplex.promote_parameter('steps', promoted_name='steps')

warmup = OpenMMwarmupNVTCube('warmup')
warmup.promote_parameter('picosec', promoted_name='warm_psec', default=10.0)

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='warmup')

job.add_cubes(ifs, idtag, smirff, complex_setup, PLmask,
              minComplex, warmup, ofs)
ifs.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)
smirff.success.connect(complex_setup.intake)
complex_setup.success.connect(PLmask.intake)
PLmask.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect( ofs.intake)

if __name__ == "__main__":
    job.run()
