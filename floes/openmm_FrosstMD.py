from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube
from OpenMMCubes.cubesPrepMD import OpenMMmakePLmaskCube, OpenMMminimizeCube, OpenMMwarmupNVTCube, OpenMMequilCube

job = WorkFloe("FrosstMD")

job.description = """
**Set up an OpenMM complex then run MD a la Merck Frosst**

This floe will do the following in each cube:
  (1) ifs: Read in the ligand file (toluene.pdb),
  (2) idtag: Add an idtag from the molecule's title or use a random 6 character string.
  (3a) smirff: Parameterize the molecule with the ffxml file (smirff99Frosst.ffxml)
        Generate the ParmEd Structure and attach it to the OEMol.
  (4) complex_setup: Paramterize the protein (T4-protein.pdb) and merge with the molecule Structure,
        Using PDBFixer: add missing atoms, add hydrogens given a pH, and solvate with TIP3P.
        Attach tagged data containing the <idtag>, <Structure> and <System>
  (5) Run restrained minimization, warmup, and equilibration similar to the Merck Frosst protocol.
  (6) Run production MD on the equilibrated complex.
  (5) ofs: Write out the OEMOl of the complex to a <idtag>-complex.oeb.gz

Ex. `python floes/openmm_FrosstMD.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb`

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
ofs: Outputs a <idtag>-complex.oeb.gz file containing the
OpenMM System, State, and ParmEd Structure of the protein:ligand complex,
packaged with the OEMol.
"""

job.classification = [['Complex Setup', 'FrosstMD']]
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

equil1 = OpenMMequilCube('equil1')
equil1.promote_parameter('picosec', promoted_name='equil1_psec', default=10.0,
      description='Length of MD run in picoseconds')
equil1.promote_parameter('restraintType', promoted_name='equil1_restrType',
      default='NonSolventNonH', description='Which kinds of atoms get xyz restraints')
equil1.promote_parameter('restraintWt', promoted_name='equil1_restrWeight',
      default=2.0, description='Restraint weight in kcal/mol/ang^2 for xyz atom restraints')
equil1.promote_parameter('label', promoted_name='equil1_label',
      default='_equil1', description='label to add to filenames from this cube')

equil2 = OpenMMequilCube('equil2')
equil2.promote_parameter('picosec', promoted_name='equil2_psec', default=20.0,
      description='Length of MD run in picoseconds')
equil2.promote_parameter('restraintType', promoted_name='equil2_restrType',
      default='NonSolventNonH', description='Which kinds of atoms get xyz restraints')
equil2.promote_parameter('restraintWt', promoted_name='equil2_restrWeight',
      default=0.5, description='Restraint weight in kcal/mol/ang^2 for xyz atom restraints')
equil2.promote_parameter('label', promoted_name='equil2_label',
      default='_equil2', description='label to add to filenames from this cube')

equil3 = OpenMMequilCube('equil3')
equil3.promote_parameter('picosec', promoted_name='equil3_psec', default=60.0,
      description='Length of MD run in picoseconds')
equil3.promote_parameter('snapFreq', promoted_name='snapFreq', default=2.0,
      description='frequency (in picoseconds) for taking snapshots')
equil3.promote_parameter('restraintType', promoted_name='equil3_restrType',
      default='CAlphaLigandNonH', description='Which kinds of atoms get xyz restraints')
equil3.promote_parameter('restraintWt', promoted_name='equil3_restrWeight',
      default=0.1, description='Restraint weight in kcal/mol/ang^2 for xyz atom restraints')
equil3.promote_parameter('label', promoted_name='equil3_label',
      default='_equil3', description='label to add to filenames from this cube')

prod = OpenMMequilCube('prod')
prod.promote_parameter('picosec', promoted_name='picosec', default=2000.0,
      description='Length of MD run in picoseconds')
prod.promote_parameter('snapFreq', promoted_name='snapFreq', default=2.0,
      description='frequency (in picoseconds) for taking snapshots')
prod.promote_parameter('restraintType', promoted_name='restraintType',
      default='None', description='Which kinds of atoms get xyz restraints')
prod.promote_parameter('label', promoted_name='label',
      default='_prod', description='label to add to filenames from this cube')

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='prod')

job.add_cubes(ifs, idtag, smirff, complex_setup, PLmask,
              minComplex, warmup, equil1, equil2, equil3, prod, ofs)
ifs.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)
smirff.success.connect(complex_setup.intake)
complex_setup.success.connect(PLmask.intake)
PLmask.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect( equil1.intake)
equil1.success.connect( equil2.intake)
equil2.success.connect( equil3.intake)
equil3.success.connect( prod.intake)
prod.success.connect( ofs.intake)

if __name__ == "__main__":
    job.run()
