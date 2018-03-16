from __future__ import unicode_literals
from floe.api import WorkFloe

from cuberecord import DataSetWriterCube
from OpenMMCubes.cubes import OpenMMminimizeSetCube, OpenMMnvtSetCube, OpenMMnptSetCube
from ComplexPrepCubes.cubes import HydrationSetCube, ComplexSetPrepCube, ForceFieldSetCube, SolvationSetCube
from ComplexPrepCubes.port import ProteinSetReaderCube
from LigPrepCubes.ports import LigandSetReaderCube, DataSetWriterCubeStripCustom
from LigPrepCubes.cubes import LigandSetChargeCube

job = WorkFloe('Merk Frosst MD Protocol')

job.description = """
Set up an OpenMM complex then minimize, warm up and equilibrate a system by using three equilibration stages

Ex: python floes/openmm_MDprep.py --ligands ligands.oeb --protein protein.oeb --ofs-data_out prep.oeb

Parameters:
-----------
ligands (file): oeb file of ligand posed in the protein active site.
protein (file): oeb file of the protein structure, assumed to be pre-prepared

Optionals:
-----------

Outputs:
--------
ofs: Outputs a ready system to MD production run
"""

job.classification = [['Complex Setup', 'FrosstMD']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = LigandSetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")


chargelig = LigandSetChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = ProteinSetReaderCube("ProteinReader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                        description="Protein file name")
iprot.promote_parameter("protein_prefix", promoted_name="protein_prefix",
                        default='prot', description="Protein prefix")

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexSetPrepCube("Complex")

# The solvation cube is used to solvate the system and define the ionic strength of the solution
# solvate = HydrationCube("Hydration")

solvate = SolvationSetCube("Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('close_solvent', promoted_name='close_solvent', default=True,
                          description='The solvent molecules will be placed very close to the solute')
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')

# Force Field Application
ff = ForceFieldSetCube("ForceField")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='amber99sbildn.xml')
ff.promote_parameter('solvent_forcefield', promoted_name='solvent_ff', default='tip3p.xml')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='GAFF2')
ff.promote_parameter('other_forcefield', promoted_name='other_ff', default='GAFF2')

# Output the prepared systems
# complex_prep_ofs = DataSetWriterCubeStripCustom('complex_prep_ofs', title='ComplexSetUpOut')
complex_prep_ofs = DataSetWriterCube('complex_prep_ofs', title='ComplexSetUpOut')
complex_prep_ofs.set_parameters(data_out=iprot.promoted_parameters['protein_prefix']['default']+'_SetUp.oeb.gz')

# Minimization
minComplex = OpenMMminimizeSetCube('minComplex', title='Minimize')
minComplex.promote_parameter('restraints', promoted_name='m_restraints', default="noh (ligand or protein)",
                             description='Select mask to apply restarints')
minComplex.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                             description='Restraint weight')
minComplex.promote_parameter('steps', promoted_name='steps', default=1000)
minComplex.promote_parameter('center', promoted_name='center', default=True)

# Output the minimized systems
# minimization_ofs = DataSetWriterCubeStripCustom('minimization_ofs', title='MinimizationOut')
minimization_ofs = DataSetWriterCube('minimization_ofs', title='MinimizationOut')

minimization_ofs.set_parameters(data_out=iprot.promoted_parameters['protein_prefix']['default']+'_Minimization.oeb.gz')

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = OpenMMnvtSetCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_psec', default=10.0,
                         description='Length of MD run in picoseconds')
warmup.promote_parameter('restraints', promoted_name='w_restraints', default="noh (ligand or protein)",
                         description='Select mask to apply restarints')
warmup.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0, description='Restraint weight')
warmup.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ps')
warmup.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=1.0,
                         description='Reporter saving intervalin ps')
warmup.promote_parameter('outfname', promoted_name='w_outfname', default='warmup',
                         description='Equilibration suffix name')

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1 = OpenMMnptSetCube('equil1', title='equil1')
equil1.promote_parameter('time', promoted_name='eq1_psec', default=10.0,
                         description='Length of MD run in picoseconds')
equil1.promote_parameter('restraints', promoted_name='eq1_restraints', default="noh (ligand or protein)",
                         description='Select mask to apply restarints')
equil1.promote_parameter('restraintWt', promoted_name='eq1_restraintWt', default=2.0, description='Restraint weight')
equil1.promote_parameter('trajectory_interval', promoted_name='eq1_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ps')
equil1.promote_parameter('reporter_interval', promoted_name='eq1_reporter_interval', default=1.0,
                         description='Reporter saving interval in ps')
equil1.promote_parameter('outfname', promoted_name='eq1_outfname', default='equil1',
                         description='Equilibration suffix name')

# NPT Equilibration stage 2
equil2 = OpenMMnptSetCube('equil2', title='equil2')
equil2.promote_parameter('time', promoted_name='eq2_psec', default=20.0,
                         description='Length of MD run in picoseconds')
equil2.promote_parameter('restraints', promoted_name='eq2_restraints', default="noh (ligand or protein)",
                         description='Select mask to apply restarints')
equil2.promote_parameter('restraintWt', promoted_name='eq2_restraintWt', default=0.5,
                         description='Restraint weight')
equil2.promote_parameter('trajectory_interval', promoted_name='eq2_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ps')
equil2.promote_parameter('reporter_interval', promoted_name='eq2_reporter_interval', default=1.0,
                         description='Reporter saving interval in ps')
equil2.promote_parameter('outfname', promoted_name='eq2_outfname', default='equil2',
                         description='Equilibration suffix name')

# NPT Equilibration stage 3
equil3 = OpenMMnptSetCube('equil3', title='equil3')
equil3.promote_parameter('time', promoted_name='eq3_psec', default=60.0,
                         description='Length of MD run in picoseconds')
equil3.promote_parameter('restraints', promoted_name='eq3_restraints', default="ca_protein or (noh ligand)",
                         description='Select mask to apply restarints')
equil3.promote_parameter('restraintWt', promoted_name='eq3_restraintWt', default=0.1,
                         description='Restraint weight')
equil3.promote_parameter('trajectory_interval', promoted_name='eq3_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ps')
equil3.promote_parameter('reporter_interval', promoted_name='eq3_reporter_interval', default=1.0,
                         description='Reporter saving interval in ps')
equil3.promote_parameter('outfname', promoted_name='eq3_outfname', default='equil3',
                         description='Equilibration suffix name')

# Output the equilibrated systems
# equilibration_ofs = DataSetWriterCubeStripCustom("equilibration_ofs", title='EquilibrationOut')
equilibration_ofs = DataSetWriterCube("equilibration_ofs", title='EquilibrationOut')
equilibration_ofs.set_parameters(data_out=iprot.promoted_parameters['protein_prefix']['default']+'_Equilibration.oeb.gz')

prod = OpenMMnptSetCube("Production")
prod.promote_parameter('time', promoted_name='prod_psec', default=2000.0,
                       description='Length of MD run in picoseconds')

prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval', default=1.0,
                       description='Trajectory saving interval in ps')
prod.promote_parameter('reporter_interval', promoted_name='prod_reporter_interval', default=1.0,
                       description='Reporter saving interval is ps')
prod.promote_parameter('outfname', promoted_name='prod_outfname', default='prod',
                       description='Equilibration suffix name')

# ofs = DataSetWriterCubeStripCustom('ofs', title='OFS-Success')
ofs = DataSetWriterCube('ofs', title='OFS-Success')


fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iprot, iligs, chargelig, complx, solvate, ff, complex_prep_ofs,
              minComplex, minimization_ofs, warmup, equil1, equil2, equil3, equilibration_ofs, prod, ofs, fail)

iprot.success.connect(complx.protein_port)
iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
complx.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
ff.success.connect(complex_prep_ofs.intake)
minComplex.success.connect(warmup.intake)
minComplex.success.connect(minimization_ofs.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(prod.intake)
equil3.success.connect(equilibration_ofs.intake)
prod.success.connect(ofs.intake)
prod.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
