#!/usr/bin/env python

# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from os import path

from floe.api import (WorkFloe,
                      ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.MDEngines.cubes import (ParallelMDMinimizeCube,
                                     ParallelMDNvtCube,
                                     ParallelMDNptCube)

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.System.cubes import (ParallelSolvationCube,
                                  MDComponentCube)

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.LigPrep.cubes import (ParallelLigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelRecordSizeCheck)

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ParallelTrajToOEMolCube,
                                                      ParallelTrajInteractionEnergyCube,
                                                      ParallelTrajPBSACube,
                                                      ConformerGatheringData,
                                                      ParallelConfTrajsToLigTraj,
                                                      ParallelConcatenateTrajMMPBSACube)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelMakeClusterTrajOEMols,
                                                       ParallelMDTrajAnalysisClusterReport,
                                                       ParallelClusterPopAnalysis,
                                                       ParallelTrajAnalysisReportDataset,
                                                       MDFloeReportCube)

job = WorkFloe('Short Trajectory MD with Analysis',
               title='Short Trajectory MD with Analysis')

job.description = open(path.join(path.dirname(__file__), 'ShortTrajMDWithAnalysis_desc.rst'), 'r').read()

job.classification = [['Specialized MD']]
job.uuid = "c831d03e-c0cb-48b0-aa02-f848da8fd1a6"
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input Dataset", description="Ligand Dataset")

ligset = LigandSetting("LigandSetting", title="Ligand Setting")
ligset.promote_parameter('max_md_runs', promoted_name='max_md_runs',
                         default=500,
                         description='The maximum allowed number of md runs')
ligset.set_parameters(lig_res_name='LIG')

chargelig = ParallelLigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)

ligid = IDSettingCube("Ligand Ids")
job.add_cube(ligid)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = DatasetReaderCube("ProteinReader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input Dataset',
                        description="Protein Dataset")

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex", title="Complex Preparation")

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = ParallelSolvationCube("Solvation", title="Solvation")

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.0')

# Protein Setting
mdcomp = MDComponentCube("MD Components", title="MD Components")
mdcomp.promote_parameter("flask_title", promoted_name="flask_title", default="")


# Production run
prod = ParallelMDNptCube("Production", title="Production")
prod.promote_parameter('time', promoted_name='prod_ns', default=2.0,
                       description='Length of MD run in nanoseconds')
prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval', default=0.004,
                       description='Trajectory saving interval in ns')
prod.promote_parameter('hmr', promoted_name="HMR", title='Use Hydrogen Mass Repartitioning', default=True,
                       description='Give hydrogens more mass to speed up the MD')
prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                       description='Select the MD Engine')
prod.set_parameters(reporter_interval=0.004)
prod.set_parameters(suffix='prod')


# Minimization
minComplex = ParallelMDMinimizeCube('minComplex', title='Minimization')
minComplex.modify_parameter(minComplex.restraints, promoted=False, default="noh (ligand or protein)")
minComplex.modify_parameter(minComplex.restraintWt, promoted=False, default=5.0)
minComplex.modify_parameter(minComplex.steps, promoted=False, default=0)
minComplex.set_parameters(center=True)
minComplex.set_parameters(save_md_stage=True)
minComplex.set_parameters(hmr=False)
minComplex.promote_parameter("md_engine", promoted_name="md_engine")

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = ParallelMDNvtCube('warmup', title='Warm Up')
warmup.set_parameters(time=0.01)
warmup.modify_parameter(warmup.restraints, promoted=False, default="noh (ligand or protein)")
warmup.modify_parameter(warmup.restraintWt, promoted=False, default=2.0)
warmup.set_parameters(trajectory_interval=0.0)
warmup.set_parameters(reporter_interval=0.001)
warmup.set_parameters(suffix='warmup')
warmup.set_parameters(hmr=False)
warmup.set_parameters(save_md_stage=True)
warmup.promote_parameter("md_engine", promoted_name="md_engine")


# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1 = ParallelMDNptCube('equil1', title='Equilibration I')
equil1.set_parameters(time=0.01)
equil1.promote_parameter("hmr", promoted_name="HMR", default=True)
equil1.modify_parameter(equil1.restraints, promoted=False, default="noh (ligand or protein)")
equil1.modify_parameter(equil1.restraintWt, promoted=False, default=1.0)
equil1.set_parameters(trajectory_interval=0.0)
equil1.set_parameters(reporter_interval=0.001)
equil1.set_parameters(suffix='equil1')
equil1.promote_parameter("md_engine", promoted_name="md_engine")


# NPT Equilibration stage 2
equil2 = ParallelMDNptCube('equil2', title='Equilibration II')
equil2.set_parameters(time=0.02)
equil2.promote_parameter("hmr", promoted_name="HMR", default=True)
equil2.modify_parameter(equil2.restraints, promoted=False, default="noh (ligand or protein)")
equil2.modify_parameter(equil2.restraintWt, promoted=False, default=0.5)
equil2.set_parameters(trajectory_interval=0.0)
equil2.set_parameters(reporter_interval=0.001)
equil2.set_parameters(suffix='equil2')
equil2.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 3
equil3 = ParallelMDNptCube('equil3', title='Equilibration III')
equil3.modify_parameter(equil3.time, promoted=False, default=0.1)
equil3.promote_parameter("hmr", promoted_name="HMR")
equil3.modify_parameter(equil3.restraints, promoted=False, default="noh (ligand or protein)")
equil3.modify_parameter(equil3.restraintWt, promoted=False, default=0.2)
equil3.set_parameters(trajectory_interval=0.0)
equil3.set_parameters(reporter_interval=0.002)
equil3.set_parameters(suffix='equil3')
equil3.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 4
equil4 = ParallelMDNptCube('equil4', title='Equilibration IV')
equil4.modify_parameter(equil4.time, promoted=False, default=0.1)
equil4.promote_parameter("hmr", promoted_name="HMR", default=True)
equil4.modify_parameter(equil4.restraints, promoted=False, default="ca_protein or (noh ligand)")
equil4.modify_parameter(equil4.restraintWt, promoted=False, default=0.1)
equil4.set_parameters(trajectory_interval=0.0)
equil4.set_parameters(reporter_interval=0.002)
equil4.set_parameters(suffix='equil4')
equil4.promote_parameter("md_engine", promoted_name="md_engine")

md_group = ParallelCubeGroup(cubes=[minComplex, warmup, equil1, equil2, equil3, equil4, prod])
job.add_group(md_group)

trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube", title="Trajectory To OEMols")
IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube", title="MM Energies")
PBSACube = ParallelTrajPBSACube("TrajPBSACube", title="PBSA Energies")

trajproc_group = ParallelCubeGroup(cubes=[trajCube, IntECube, PBSACube])
job.add_group(trajproc_group)

confGather = ConformerGatheringData("Gathering Conformer Records", title="Gathering Conformer Records")
catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj", title="Combine Pose Trajectories")
catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube', title="Concatenate MMPBSA Energies")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube", title="Clustering")
clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis', title="Clustering Analysis")
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols', title="Per-Cluster Analysis")
prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset', title="Analysis Report")
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport", title="Relevant Output Extraction")

analysis_group = ParallelCubeGroup(cubes=[catLigTraj, catLigMMPBSA, clusCube, clusPop,
                                          clusOEMols, prepDataset, report_gen])
job.add_group(analysis_group)

report = MDFloeReportCube("report", title="Floe Report")

# This cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success", title="Record Check Success")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(iligs, ligset, iprot, mdcomp, chargelig, complx,
              solvate, coll_open, ff,
              minComplex, warmup, equil1, equil2, equil3, equil4, prod,
              trajCube, IntECube, PBSACube, confGather,
              catLigTraj, catLigMMPBSA, clusCube, clusPop, clusOEMols,
              prepDataset, report_gen, report,
              coll_close, check_rec, ofs, fail)

# Success Connections
iligs.success.connect(ligset.intake)
ligset.success.connect(chargelig.intake)
chargelig.success.connect(ligid.intake)
ligid.success.connect(complx.intake)
iprot.success.connect(mdcomp.intake)
mdcomp.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)
solvate.success.connect(coll_open.intake)
coll_open.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(equil4.intake)
equil4.success.connect(prod.intake)
prod.success.connect(trajCube.intake)
trajCube.success.connect(IntECube.intake)
IntECube.success.connect(PBSACube.intake)
PBSACube.success.connect(confGather.intake)
confGather.success.connect(catLigTraj.intake)
catLigTraj.success.connect(catLigMMPBSA.intake)
catLigMMPBSA.success.connect(clusCube.intake)
clusCube.success.connect(clusPop.intake)
clusPop.success.connect(clusOEMols.intake)
clusOEMols.success.connect(prepDataset.intake)
prepDataset.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail Connections
ligset.failure.connect(check_rec.fail_in)
chargelig.failure.connect(check_rec.fail_in)
ligid.failure.connect(check_rec.fail_in)
mdcomp.failure.connect(check_rec.fail_in)
complx.failure.connect(check_rec.fail_in)
solvate.failure.connect(check_rec.fail_in)
coll_open.failure.connect(check_rec.fail_in)
ff.failure.connect(check_rec.fail_in)
minComplex.failure.connect(check_rec.fail_in)
warmup.failure.connect(check_rec.fail_in)
equil1.failure.connect(check_rec.fail_in)
equil2.failure.connect(check_rec.fail_in)
equil3.failure.connect(check_rec.fail_in)
equil4.failure.connect(check_rec.fail_in)
prod.failure.connect(check_rec.fail_in)
trajCube.failure.connect(check_rec.fail_in)
IntECube.failure.connect(check_rec.fail_in)
PBSACube.failure.connect(check_rec.fail_in)
confGather.failure.connect(check_rec.fail_in)
catLigTraj.failure.connect(check_rec.fail_in)
catLigMMPBSA.failure.connect(check_rec.fail_in)
clusCube.failure.connect(check_rec.fail_in)
clusPop.failure.connect(check_rec.fail_in)
clusOEMols.failure.connect(check_rec.fail_in)
prepDataset.failure.connect(check_rec.fail_in)
report_gen.failure.connect(check_rec.fail_in)
report.failure.connect(check_rec.fail_in)
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
