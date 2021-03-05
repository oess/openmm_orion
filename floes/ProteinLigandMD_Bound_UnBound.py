from floe.api import (WorkFloe,
                      ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.LigPrep.cubes import (ParallelLigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelSolvationCube)

from MDOrion.MDEngines.cubes import (ParallelMDMinimizeCube,
                                     ParallelMDNvtCube,
                                     ParallelMDNptCube)

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.System.cubes import MDComponentCube

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.FEC.RFEC.cubes import (BoundUnboundSwitchCube,
                                    RBFECMapping)

job = WorkFloe("MD of Ligand Bound and Unbound to Protein",
               title="MD of Ligand Bound and Unbound to Protein")

job.description = """
TBD
"""

job.classification = [['MD', 'Protein-Ligand', 'FEC']]
job.uuid = "537f64c5-0d84-4537-ad74-c55037304e07"
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input Dataset", description="Ligand Dataset")

ligset = LigandSetting("Ligand Setting", title="Ligand Setting")
ligset.set_parameters(lig_res_name='LIG')

chargelig = ParallelLigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)

ligid = IDSettingCube("Ligand Ids")

md_lig_components = MDComponentCube("MD Ligand Components", title="MD Ligand Components")
md_lig_components.set_parameters(flask_title="")
md_lig_components.set_parameters(multiple_flasks=True)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = DatasetReaderCube("ProteinReader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input Dataset',
                        description="Protein Dataset")

md_prot_components = MDComponentCube("MD Protein Components", title="MD Protein Components")
md_prot_components.promote_parameter("flask_title", promoted_name="flask_title",
                                     description='Prefix name used to identity the Protein', default='')
md_prot_components.set_parameters(multiple_flasks=False)

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex", title="Complex Preparation")

solvate = ParallelSolvationCube("Solvation", title="Solvation")
solvate.set_parameters(density=1.03)
solvate.set_parameters(salt_concentration=50.0)
solvate.modify_parameter(solvate.close_solvent, promoted=False, default=False)

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)
coll_open.set_parameters(write_new_collection='MD_OPLMD')

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.0')

# Switching Bound and Unbound runs
switch = BoundUnboundSwitchCube("Bound/Unbound In Switch", title='Bound/Unbound In Switch')

# Run the equilibrium Simulations ot the Unbound-States
prod_uns = ParallelMDNptCube("Production Unbound States", title="Production Unbound States")
prod_uns.promote_parameter('time', promoted_name='prod_us_ns', default=6.0,
                           description='Length of Production MD run in nanoseconds')
# prod_uns.promote_parameter('trajectory_frames', promoted_name='prod_trajectory_us_frames', default=1500,
#                            description='Total number of trajectory frames used in the NES calculation')
prod_uns.modify_parameter(prod_uns.trajectory_frames, promoted=False, default=1500)
prod_uns.promote_parameter('hmr', promoted_name="hmr_us", title='Use Hydrogen Mass Repartitioning '
                                                                'in the Unbound simulation', default=True,
                           description='Give hydrogens more mass to speed up the MD')
prod_uns.set_parameters(reporter_interval=0.002)
prod_uns.set_parameters(suffix='prod_unb')
prod_uns.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM', description='Select the MD Engine')

# Minimization of the Unbound-States
minimize_uns = ParallelMDMinimizeCube("Minimize Unbound States", title="Minimization Unbound States")
minimize_uns.set_parameters(restraints='noh ligand')
minimize_uns.promote_parameter("md_engine", promoted_name="md_engine")
minimize_uns.set_parameters(steps=2000)
minimize_uns.set_parameters(restraintWt=5.0)
minimize_uns.set_parameters(center=True)
minimize_uns.set_parameters(hmr=False)

# NVT Warm-up of the Unbound-States
warmup_uns = ParallelMDNvtCube('Warmup Unbound States', title='Warmup Unbound States')
warmup_uns.set_parameters(time=1)
warmup_uns.set_parameters(restraints="noh ligand")
warmup_uns.promote_parameter("md_engine", promoted_name="md_engine")
warmup_uns.set_parameters(restraintWt=2.0)
warmup_uns.set_parameters(trajectory_interval=0.0)
warmup_uns.set_parameters(reporter_interval=0.002)
warmup_uns.set_parameters(suffix='warmup_un')
warmup_uns.set_parameters(hmr=False)

# NPT Equilibration stage of the Unbound-States
equil_uns = ParallelMDNptCube('Equilibration Unbound States', title='Equilibration Unbound States')
equil_uns.set_parameters(time=1)
equil_uns.promote_parameter("hmr", promoted_name="hmr_us", default=True)
equil_uns.set_parameters(restraints="noh ligand")
equil_uns.promote_parameter("md_engine", promoted_name="md_engine")
equil_uns.set_parameters(restraintWt=0.1)
equil_uns.set_parameters(trajectory_interval=0.0)
equil_uns.set_parameters(reporter_interval=0.002)
equil_uns.set_parameters(suffix='equil_un')

md_group_uns = ParallelCubeGroup(cubes=[minimize_uns, warmup_uns, equil_uns, prod_uns])
job.add_group(md_group_uns)

# Run the equilibrium Simulations ot the Bound-States
prod_bns = ParallelMDNptCube("Production Bound States", title="Production Bound States")
prod_bns.promote_parameter('time', promoted_name='prod_us_ns', default=6.0)
# prod_bns.promote_parameter('trajectory_frames', promoted_name='prod_trajectory_us_frames', default=1500)
prod_bns.modify_parameter(prod_bns.trajectory_frames, promoted=False, default=1500)
prod_bns.promote_parameter('hmr', promoted_name="hmr_bs", title='Use Hydrogen Mass Repartitioning '
                                                                'in the Bound simulation', default=True,
                           description='Give hydrogens more mass to speed up the MD')
prod_bns.modify_parameter(prod_bns.trajectory_frames, promoted=False, default=1500)
prod_bns.promote_parameter("md_engine", promoted_name="md_engine")
prod_bns.set_parameters(reporter_interval=0.002)
prod_bns.set_parameters(suffix='prod_bn')

# Minimization of the Bound-States
minimize_bns = ParallelMDMinimizeCube("Minimize Bound States", title="Minimization Bound States")
minimize_bns.modify_parameter(minimize_bns.restraints, promoted=False, default="noh (ligand or protein)")
minimize_bns.modify_parameter(minimize_bns.restraintWt, promoted=False, default=5.0)
minimize_bns.promote_parameter("md_engine", promoted_name="md_engine")
minimize_bns.set_parameters(steps=2000)
minimize_bns.set_parameters(center=True)
minimize_bns.set_parameters(save_md_stage=True)
minimize_bns.set_parameters(hmr=False)

# NVT Warm-up of the Unbound-States
warmup_bns = ParallelMDNvtCube('Warmup Bound States', title='Warmup Bound States')
warmup_bns.set_parameters(time=0.01)
warmup_bns.modify_parameter(warmup_bns.restraints, promoted=False, default="noh (ligand or protein)")
warmup_bns.modify_parameter(warmup_bns.restraintWt, promoted=False, default=2.0)
warmup_bns.promote_parameter("md_engine", promoted_name="md_engine")
warmup_bns.set_parameters(trajectory_interval=0.0)
warmup_bns.set_parameters(reporter_interval=0.001)
warmup_bns.set_parameters(suffix='warmup_bn')
warmup_bns.set_parameters(hmr=False)
warmup_bns.set_parameters(save_md_stage=True)

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Bound Equilibration stage 1
equil1_bns = ParallelMDNptCube('equil1_bns', title='Equilibration I Bound States')
equil1_bns.set_parameters(time=0.01)
equil1_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil1_bns.modify_parameter(equil1_bns.restraints, promoted=False, default="noh (ligand or protein)")
equil1_bns.modify_parameter(equil1_bns.restraintWt, promoted=False, default=1.0)
equil1_bns.promote_parameter("md_engine", promoted_name="md_engine")
equil1_bns.set_parameters(trajectory_interval=0.0)
equil1_bns.set_parameters(reporter_interval=0.001)
equil1_bns.set_parameters(suffix='equil1_bn')

# NPT Bound Equilibration stage 2
equil2_bns = ParallelMDNptCube('equil2_bs', title='Equilibration II Bound States')
equil2_bns.set_parameters(time=0.02)
equil2_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil2_bns.modify_parameter(equil2_bns.restraints, promoted=False, default="noh (ligand or protein)")
equil2_bns.modify_parameter(equil2_bns.restraintWt, promoted=False, default=0.5)
equil2_bns.promote_parameter("md_engine", promoted_name="md_engine")
equil2_bns.set_parameters(trajectory_interval=0.0)
equil2_bns.set_parameters(reporter_interval=0.001)
equil2_bns.set_parameters(suffix='equil2_bs')

# NPT Bound Equilibration stage 3
equil3_bns = ParallelMDNptCube('equil3_bs', title='Equilibration III Bound States')
equil3_bns.set_parameters(time=0.1)
equil3_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil3_bns.modify_parameter(equil3_bns.restraints, promoted=False, default="ca_protein or (noh ligand)")
equil3_bns.modify_parameter(equil3_bns.restraintWt, promoted=False, default=0.2)
equil3_bns.promote_parameter("md_engine", promoted_name="md_engine")
equil3_bns.set_parameters(trajectory_interval=0.0)
equil3_bns.set_parameters(reporter_interval=0.002)
equil3_bns.set_parameters(suffix='equil3_bn')

# NPT Equilibration stage 4
equil4_bns = ParallelMDNptCube('equil4_bs', title='Equilibration IV Bound States')
equil4_bns.modify_parameter(equil4_bns.time, promoted=False, default=0.1)
equil4_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil4_bns.modify_parameter(equil4_bns.restraints, promoted=False, default="ca_protein or (noh ligand)")
equil4_bns.modify_parameter(equil4_bns.restraintWt, promoted=False, default=0.1)
equil4_bns.promote_parameter("md_engine", promoted_name="md_engine")
equil4_bns.set_parameters(trajectory_interval=0.0)
equil4_bns.set_parameters(reporter_interval=0.002)
equil4_bns.set_parameters(suffix='equil4_bn')

md_group_bs = ParallelCubeGroup(cubes=[minimize_bns, warmup_bns, equil1_bns, equil2_bns, equil3_bns, equil4_bns, prod_bns])
job.add_group(md_group_bs)


# This cube is necessary for the correct working of collections and shards
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

# Switching Bound and Unbound runs
switch_out = BoundUnboundSwitchCube("Bound/Unbound Out Switch", title='Bound/Unbound Out Switch')

ofs_unbound = DatasetWriterCube('ofs_unbound', title='Unbound Out')
ofs_unbound.promote_parameter("data_out", promoted_name="out_unbound", title="Unbound Dataset Out", description="Unbound Dataset Out")

ofs_bound = DatasetWriterCube('ofs_bound', title='Bound Out')
ofs_bound.promote_parameter("data_out", promoted_name="out_bound", title="Bound Dataset Out", description="Bound Dataset Out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="Dataset Failures out")

job.add_cubes(iligs, ligset, chargelig, ligid, md_lig_components, coll_open,
              iprot, md_prot_components, complx, solvate, ff, switch,
              minimize_bns, warmup_bns, equil1_bns,
              equil2_bns, equil3_bns, equil4_bns, prod_bns,
              minimize_uns, warmup_uns, equil_uns, prod_uns,
              coll_close, switch_out, fail, ofs_unbound, ofs_bound)

# Ligand Setting
iligs.success.connect(ligset.intake)
ligset.success.connect(chargelig.intake)
chargelig.success.connect(ligid.intake)
ligid.success.connect(md_lig_components.intake)
md_lig_components.success.connect(coll_open.intake)
coll_open.success.connect(complx.intake)
coll_open.success.connect(solvate.intake)

# Protein Setting
iprot.success.connect(md_prot_components.intake)
md_prot_components.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)

# Flask Setting
solvate.success.connect(ff.intake)
ff.success.connect(switch.intake)

# Unbound MD Run
switch.success.connect(minimize_uns.intake)
minimize_uns.success.connect(warmup_uns.intake)
warmup_uns.success.connect(equil_uns.intake)
equil_uns.success.connect(prod_uns.intake)
prod_uns.success.connect(coll_close.intake)

# Bound MD run
switch.bound_port.connect(minimize_bns.intake)
minimize_bns.success.connect(warmup_bns.intake)
warmup_bns.success.connect(equil1_bns.intake)
equil1_bns.success.connect(equil2_bns.intake)
equil2_bns.success.connect(equil3_bns.intake)
equil3_bns.success.connect(equil4_bns.intake)
equil4_bns.success.connect(prod_bns.intake)
prod_bns.success.connect(coll_close.intake)

coll_close.success.connect(switch_out.intake)
switch_out.success.connect(ofs_unbound.intake)
switch_out.bound_port.connect(ofs_bound.intake)

# Fail port connections
ligset.failure.connect(fail.intake)
chargelig.failure.connect(fail.intake)
ligid.failure.connect(fail.intake)
md_lig_components.failure.connect(fail.intake)
md_prot_components.failure.connect(fail.intake)
complx.failure.connect(fail.intake)
coll_open.failure.connect(fail.intake)
solvate.failure.connect(fail.intake)
ff.failure.connect(fail.intake)
switch.failure.connect(fail.intake)

minimize_uns.failure.connect(fail.intake)
warmup_uns.failure.connect(fail.intake)
equil_uns.failure.connect(fail.intake)
prod_uns.failure.connect(fail.intake)

minimize_bns.failure.connect(fail.intake)
warmup_bns.failure.connect(fail.intake)
equil1_bns.failure.connect(fail.intake)
equil2_bns.failure.connect(fail.intake)
equil3_bns.failure.connect(fail.intake)
equil4_bns.failure.connect(fail.intake)
prod_bns.failure.connect(fail.intake)

coll_close.failure.connect(fail.intake)
switch_out.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()


