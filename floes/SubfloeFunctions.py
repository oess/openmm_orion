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
from MDOrion.System.cubes import (ParallelRecordSizeCheck)
from MDOrion.System.cubes import CollectionSetting

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

from MDOrion.TrjAnalysis.cubes_hintAnalysis import (ParallelComparePoseBintsToTrajBints)


def setup_MD_startup(input_floe, input_cube, output_cube, fail_cube, options):
    # Production run
    prod = ParallelMDNptCube("Production", title="Production")
    prod.promote_parameter('time', promoted_name='prod_ns',
                           default=options['Prod_Default_Time_ns'],
                           description='Length of MD run in nanoseconds')
    prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval',
                           default=options['Prod_Default_Traj_Intvl_ns'],
                           description='Trajectory saving interval in ns')
    prod.promote_parameter('hmr', promoted_name="HMR", title='Use Hydrogen Mass Repartitioning', default=True,
                           description='Give hydrogens more mass to speed up the MD')
    prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                           description='Select the MD Engine')
    prod.set_parameters(reporter_interval=options['Prod_Default_Traj_Intvl_ns'])
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

    # The system is equilibrated at the right pressure and temperature in several stages
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
    equil3.promote_parameter("hmr", promoted_name="HMR", default=True)
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
    input_floe.add_group(md_group)

    input_floe.add_cubes(minComplex, warmup, equil1, equil2, equil3, equil4, prod)

    # Success Connections
    input_cube.success.connect(minComplex.intake)
    minComplex.success.connect(warmup.intake)
    warmup.success.connect(equil1.intake)
    equil1.success.connect(equil2.intake)
    equil2.success.connect(equil3.intake)
    equil3.success.connect(equil4.intake)
    equil4.success.connect(prod.intake)
    prod.success.connect(output_cube.intake)
    
    # Fail Connections
    minComplex.failure.connect(fail_cube.fail_in)
    warmup.failure.connect(fail_cube.fail_in)
    equil1.failure.connect(fail_cube.fail_in)
    equil2.failure.connect(fail_cube.fail_in)
    equil3.failure.connect(fail_cube.fail_in)
    equil4.failure.connect(fail_cube.fail_in)
    prod.failure.connect(fail_cube.fail_in)

    return True


def setup_traj_analysis(input_floe, input_cube, output_cube, fail_cube):
    trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube", title="Trajectory To OEMols")
    IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube", title="MM Energies")
    PBSACube = ParallelTrajPBSACube("TrajPBSACube", title="PBSA Energies")

    trajproc_group = ParallelCubeGroup(cubes=[trajCube, IntECube, PBSACube])
    input_floe.add_group(trajproc_group)

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
    input_floe.add_group(analysis_group)

    report = MDFloeReportCube("report", title="Floe Report")

    input_floe.add_cubes(trajCube, IntECube, PBSACube, confGather,
                  catLigTraj, catLigMMPBSA, clusCube, clusPop, clusOEMols,
                  prepDataset, report_gen, report)
    
    # Success Connections
    input_cube.success.connect(trajCube.intake)
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
    report.success.connect(output_cube.intake)
    
    # Fail Connections
    trajCube.failure.connect(fail_cube.fail_in)
    IntECube.failure.connect(fail_cube.fail_in)
    PBSACube.failure.connect(fail_cube.fail_in)
    confGather.failure.connect(fail_cube.fail_in)
    catLigTraj.failure.connect(fail_cube.fail_in)
    catLigMMPBSA.failure.connect(fail_cube.fail_in)
    clusCube.failure.connect(fail_cube.fail_in)
    clusPop.failure.connect(fail_cube.fail_in)
    clusOEMols.failure.connect(fail_cube.fail_in)
    prepDataset.failure.connect(fail_cube.fail_in)
    report_gen.failure.connect(fail_cube.fail_in)
    report.failure.connect(fail_cube.fail_in)

    return True


def setup_bint(input_floe, input_cube, output_cube, fail_cube):

    trajBints = ParallelComparePoseBintsToTrajBints("TrajBintsCube", title="Trajectory Binding Interactions")

    input_floe.add_cubes(trajBints)

    input_cube.success.connect(trajBints.intake)
    trajBints.success.connect(output_cube.intake)

    trajBints.failure.connect(fail_cube.fail_in)

    return True


if __name__ == "__main__":
    job.run()
