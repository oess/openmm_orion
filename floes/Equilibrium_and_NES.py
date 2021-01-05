from floe.api import (WorkFloe,
                      ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.LigPrep.cubes import (ParallelLigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelSolvationCube,
                                  ParallelRecordSizeCheck)

from MDOrion.MDEngines.cubes import (ParallelMDMinimizeCube,
                                     ParallelMDNvtCube,
                                     ParallelMDNptCube)

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.System.cubes import MDComponentCube

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.FEC.RFEC.cubes import (BoundUnboundSwitchCube,
                                    RBFECMapping,
                                    RBFECEdgeGathering,
                                    ParallelGMXChimera,
                                    ParallelNESGMX,
                                    NESAnalysis)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import MDFloeReportCube

job = WorkFloe("Equilibrium and Non Equilibrium Switching", title="Equilibrium and Non Equilibrium Switching")

job.description = """
The Non Equilibrium Switching (NES) floe performs Relative Binding Affinity Calculations
between a set of provided ligands and the relate provided target protein. The floe
requires also a text file with the map of the edges between the different relative
free energy calculations to run. The file format of the text file is a set of lines
with the syntax:

ligA_name >> ligB_name

where ligA_name and ligB_name are respectively strings of the ligand names for the 
ligand in the starting state A and  the ligand name in the final state B. Because the 
edges to run are defined by ligand names it is important that all the submitted ligands 
have unique names. At the end of the calculations the NES floe will produce a floe 
report where the insight of each edge calculation is reported with different metrics 
used to compute the relative binding affinity.
"""

job.classification = [['FEC']]
job.uuid = "45776760-785d-4128-972e-1be13baddfc0"
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input Dataset", description="Ligand Dataset")

ligset = LigandSetting("Ligand Setting", title="Ligand Setting")
ligset.set_parameters(lig_res_name='LIG')

rbfec_map = RBFECMapping("Edge Mapping", title="Edge Mapping")
rbfec_map.promote_parameter("map_file", promoted_name="map")

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

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex", title="Complex Preparation")

solvate = ParallelSolvationCube("Solvation", title="Solvation")
solvate.set_parameters(density=1.03)
solvate.set_parameters(salt_concentration=50.0)
solvate.modify_parameter(solvate.close_solvent, promoted=False, default=False)

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.0')

# Switching Bound and Unbound runs
switch = BoundUnboundSwitchCube("Bound/Unbound Switch", title='Bound/Unbound Switch')

# Run the equilibrium Simulations ot the Unbound-States
prod_uns = ParallelMDNptCube("Production Unbound States", title="Production Unbound States")
prod_uns.promote_parameter('time', promoted_name='prod_us_ns', default=6.0,
                           description='Length of Unbound MD run in nanoseconds')
prod_uns.promote_parameter('trajectory_frames', promoted_name='prod_trajectory_us_frames', default=80,
                           description='Total number of trajectory frames used in the NES calculation')
prod_uns.promote_parameter('hmr', promoted_name="hmr_us", title='Use Hydrogen Mass Repartitioning '
                                                                'in the Unbound simulation', default=True,
                           description='Give hydrogens more mass to speed up the MD')
prod_uns.set_parameters(md_engine='OpenMM')
prod_uns.set_parameters(reporter_interval=0.002)
prod_uns.set_parameters(suffix='prod_unb')

# Minimization of the Unbound-States
minimize_uns = ParallelMDMinimizeCube("Minimize Unbound States", title="Minimization Unbound States")
minimize_uns.set_parameters(restraints='noh ligand')
minimize_uns.set_parameters(md_engine='OpenMM')
minimize_uns.set_parameters(steps=2000)
minimize_uns.set_parameters(restraintWt=5.0)
minimize_uns.set_parameters(center=True)
minimize_uns.set_parameters(hmr=False)

# NVT Warm-up of the Unbound-States
warmup_uns = ParallelMDNvtCube('Warmup Unbound States', title='Warmup Unbound States')
warmup_uns.set_parameters(time=1)
warmup_uns.set_parameters(restraints="noh ligand")
warmup_uns.set_parameters(md_engine='OpenMM')
warmup_uns.set_parameters(restraintWt=2.0)
warmup_uns.set_parameters(trajectory_interval=0.0)
warmup_uns.set_parameters(reporter_interval=0.002)
warmup_uns.set_parameters(suffix='warmup_un')
warmup_uns.set_parameters(hmr=False)

# NPT Equilibration stage of the Unbound-States
equil_uns = ParallelMDNptCube('Equilibration Unbond States', title='Equilibration Unbond States')
equil_uns.set_parameters(time=1)
equil_uns.promote_parameter("hmr", promoted_name="hmr_us", default=True)
equil_uns.set_parameters(restraints="noh ligand")
equil_uns.set_parameters(md_engine='OpenMM')
equil_uns.set_parameters(restraintWt=0.1)
equil_uns.set_parameters(trajectory_interval=0.0)
equil_uns.set_parameters(reporter_interval=0.002)
equil_uns.set_parameters(suffix='equil_un')

md_group_uns = ParallelCubeGroup(cubes=[minimize_uns, warmup_uns, equil_uns, prod_uns])
job.add_group(md_group_uns)

# Run the equilibrium Simulations ot the Bound-States
prod_bns = ParallelMDNptCube("Production Bound States", title="Production Bound States")
prod_bns.promote_parameter('time', promoted_name='prod_us_ns', default=6.0)
prod_bns.promote_parameter('trajectory_frames', promoted_name='prod_trajectory_us_frames', default=80)
prod_bns.promote_parameter('hmr', promoted_name="hmr_bs", title='Use Hydrogen Mass Repartitioning '
                                                                'in the Bound simulation', default=True,
                           description='Give hydrogens more mass to speed up the MD')
prod_bns.set_parameters(md_engine='OpenMM')
prod_bns.set_parameters(reporter_interval=0.002)
prod_bns.set_parameters(suffix='prod_bn')

# Minimization of the Bound-States
minimize_bns = ParallelMDMinimizeCube("Minimize Bound States", title="Minimization Bound States")
minimize_bns.modify_parameter(minimize_bns.restraints, promoted=False, default="noh (ligand or protein)")
minimize_bns.modify_parameter(minimize_bns.restraintWt, promoted=False, default=5.0)
minimize_bns.set_parameters(md_engine='OpenMM')
minimize_bns.set_parameters(steps=2000)
minimize_bns.set_parameters(center=True)
minimize_bns.set_parameters(save_md_stage=True)
minimize_bns.set_parameters(hmr=False)
minimize_bns.set_parameters()

# NVT Warm-up of the Unbound-States
warmup_bns = ParallelMDNvtCube('Warmup Bound States', title='Warmup Bound States')
warmup_bns.set_parameters(time=0.01)
warmup_bns.modify_parameter(warmup_bns.restraints, promoted=False, default="noh (ligand or protein)")
warmup_bns.modify_parameter(warmup_bns.restraintWt, promoted=False, default=2.0)
warmup_bns.set_parameters(md_engine='OpenMM')
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
equil_uns.set_parameters(md_engine='OpenMM')
equil1_bns.set_parameters(trajectory_interval=0.0)
equil1_bns.set_parameters(reporter_interval=0.001)
equil1_bns.set_parameters(suffix='equil1_bn')

# NPT Bound Equilibration stage 2
equil2_bns = ParallelMDNptCube('equil2_bs', title='Equilibration II Bound States')
equil2_bns.set_parameters(time=0.02)
equil2_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil2_bns.modify_parameter(equil2_bns.restraints, promoted=False, default="noh (ligand or protein)")
equil2_bns.modify_parameter(equil2_bns.restraintWt, promoted=False, default=0.5)
equil2_bns.set_parameters(md_engine='OpenMM')
equil2_bns.set_parameters(trajectory_interval=0.0)
equil2_bns.set_parameters(reporter_interval=0.001)
equil2_bns.set_parameters(suffix='equil2_bs')

# NPT Bound Equilibration stage 3
equil3_bns = ParallelMDNptCube('equil3_bs', title='Equilibration III Bound States')
equil3_bns.set_parameters(time=0.1)
equil3_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil3_bns.modify_parameter(equil3_bns.restraints, promoted=False, default="ca_protein or (noh ligand)")
equil3_bns.modify_parameter(equil3_bns.restraintWt, promoted=False, default=0.2)
equil3_bns.set_parameters(md_engine='OpenMM')
equil3_bns.set_parameters(trajectory_interval=0.0)
equil3_bns.set_parameters(reporter_interval=0.002)
equil3_bns.set_parameters(suffix='equil3_bn')

# NPT Equilibration stage 4
equil4_bns = ParallelMDNptCube('equil4_bs', title='Equilibration IV Bound States')
equil4_bns.modify_parameter(equil4_bns.time, promoted=False, default=0.1)
equil4_bns.promote_parameter("hmr", promoted_name="hmr_bs", default=True)
equil4_bns.modify_parameter(equil4_bns.restraints, promoted=False, default="ca_protein or (noh ligand)")
equil4_bns.modify_parameter(equil4_bns.restraintWt, promoted=False, default=0.1)
equil4_bns.set_parameters(md_engine='OpenMM')
equil4_bns.set_parameters(trajectory_interval=0.0)
equil4_bns.set_parameters(reporter_interval=0.002)
equil4_bns.set_parameters(suffix='equil4_bn')

md_group_bs = ParallelCubeGroup(cubes=[minimize_bns, warmup_bns, equil1_bns, equil2_bns, equil3_bns, equil4_bns, prod_bns])
job.add_group(md_group_bs)

gathering = RBFECEdgeGathering("Gathering", title="Gathering Equilibrium Runs")
gathering.promote_parameter('map_file', promoted_name='map')

chimera = ParallelGMXChimera("GMXChimera", title="GMX Chimera")

unbound_eq_nes = ParallelNESGMX("GMXUnboundEQ", title="GMX Unbound NPT Equilibration")
unbound_eq_nes.set_parameters(time=0.01)
unbound_eq_nes.set_parameters(enable_switching=False)

unbound_nes = ParallelNESGMX("GMXUnboundNES", title="GMX Unbound NES")
unbound_nes.promote_parameter("time", promoted_name="nes_time", default=0.05)
unbound_nes.set_parameters(enable_switching=True)

md_group_unbound_nes = ParallelCubeGroup(cubes=[unbound_eq_nes, unbound_nes])
job.add_group(md_group_unbound_nes)

bound_eq_nes = ParallelNESGMX("GMXBoundEQ", title="GMX Bound NPT Equilibration")
bound_eq_nes.set_parameters(time=0.02)
bound_eq_nes.set_parameters(enable_switching=False)

bound_nes = ParallelNESGMX("GMXBoundNES", title="GMX Bound NES")
bound_nes.promote_parameter("time", promoted_name="nes_time")
bound_nes.set_parameters(enable_switching=True)

md_group_bound_nes = ParallelCubeGroup(cubes=[bound_eq_nes, bound_nes])
job.add_group(md_group_bound_nes)

# This cube is necessary for the correct working of collections and shards
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

nes_analysis = NESAnalysis("NES_Analysis")

report = MDFloeReportCube("report", title="Floe Report")
report.set_parameters(floe_report_title="NES Report")

check_rec = ParallelRecordSizeCheck("Record Check Success")

ofs = DatasetWriterCube('ofs', title='NES Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="NES Dataset Out",
                      description="NES Dataset Out")

fail = DatasetWriterCube('fail', title='NES Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="NES Failures",
                       description="NES Dataset Failures out")


# TODO DEBUG ONLY TO BE REMOVED
ofs_lig = DatasetWriterCube('ofs_unbound', title='Equilibrium Unbound Out')
ofs_lig.promote_parameter("data_out", promoted_name="out_unbound",
                          title="Equilibrium Unbound Out",
                          description="Equilibrium Unbound Out")

ofs_prot = DatasetWriterCube('ofs_bound', title='Equilibrium Bound Out')
ofs_prot.promote_parameter("data_out", promoted_name="out_bound",
                           title="Equilibrium Bound Out",
                           description="Equilibrium Bound Out")

# TODO DEBUGGING ONLY REMOVE EXTRA OFS AT THE END
job.add_cubes(iligs, ligset, rbfec_map, chargelig, ligid, md_lig_components, coll_open,
              iprot, md_prot_components, complx, solvate, ff, switch,
              minimize_uns, warmup_uns, equil_uns, prod_uns,
              minimize_bns, warmup_bns, equil1_bns,
              equil2_bns, equil3_bns, equil4_bns, prod_bns, gathering,
              chimera, unbound_eq_nes, unbound_nes,
              bound_eq_nes, bound_nes,
              coll_close, nes_analysis, report,
              check_rec, ofs, fail, ofs_lig, ofs_prot)

# Ligand Setting
iligs.success.connect(ligset.intake)
ligset.success.connect(rbfec_map.intake)
rbfec_map.success.connect(chargelig.intake)
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
prod_uns.success.connect(gathering.intake)
# TODO DEBUG ONLY
prod_uns.success.connect(ofs_lig.intake)

# Bound MD run
switch.bound_port.connect(minimize_bns.intake)
minimize_bns.success.connect(warmup_bns.intake)
warmup_bns.success.connect(equil1_bns.intake)
equil1_bns.success.connect(equil2_bns.intake)
equil2_bns.success.connect(equil3_bns.intake)
equil3_bns.success.connect(equil4_bns.intake)
equil4_bns.success.connect(prod_bns.intake)
prod_bns.success.connect(gathering.bound_port)
# TODO DEBUG ONLY
prod_bns.success.connect(ofs_prot.intake)

# Chimera NES Setting
gathering.success.connect(chimera.intake)

chimera.success.connect(unbound_eq_nes.intake)
unbound_eq_nes.success.connect(unbound_nes.intake)
unbound_nes.success.connect(coll_close.intake)

chimera.bound_port.connect(bound_eq_nes.intake)
bound_eq_nes.success.connect(bound_nes.intake)
bound_nes.success.connect(coll_close.intake)

coll_close.success.connect(nes_analysis.intake)
nes_analysis.success.connect(report.intake)
report.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail port connections
ligset.failure.connect(check_rec.fail_in)
rbfec_map.failure.connect(check_rec.fail_in)
chargelig.failure.connect(check_rec.fail_in)
ligid.failure.connect(check_rec.fail_in)
md_lig_components.failure.connect(check_rec.fail_in)
md_prot_components.failure.connect(check_rec.fail_in)
complx.failure.connect(check_rec.fail_in)
coll_open.failure.connect(check_rec.fail_in)
solvate.failure.connect(check_rec.fail_in)
ff.failure.connect(check_rec.fail_in)
switch.failure.connect(check_rec.fail_in)

minimize_uns.failure.connect(check_rec.fail_in)
warmup_uns.failure.connect(check_rec.fail_in)
equil_uns.failure.connect(check_rec.fail_in)
prod_uns.failure.connect(check_rec.fail_in)

minimize_bns.failure.connect(check_rec.fail_in)
warmup_bns.failure.connect(check_rec.fail_in)
equil1_bns.failure.connect(check_rec.fail_in)
equil2_bns.failure.connect(check_rec.fail_in)
equil3_bns.failure.connect(check_rec.fail_in)
equil4_bns.failure.connect(check_rec.fail_in)
prod_bns.failure.connect(check_rec.fail_in)

gathering.failure.connect(check_rec.fail_in)
chimera.failure.connect(check_rec.fail_in)
unbound_eq_nes.failure.connect(check_rec.fail_in)
unbound_nes.failure.connect(check_rec.fail_in)
bound_eq_nes.failure.connect(check_rec.fail_in)
bound_nes.failure.connect(check_rec.fail_in)

coll_close.failure.connect(check_rec.fail_in)
nes_analysis.failure.connect(check_rec.fail_in)
report.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()


