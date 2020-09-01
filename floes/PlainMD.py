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

from MDOrion.System.cubes import MDComponentCube

from MDOrion.System.cubes import ParallelSolvationCube

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelRecordSizeCheck)

job = WorkFloe('Solvate and Run MD',
               title='Solvate and Run MD')

job.description = open(path.join(path.dirname(__file__), 'PlainMD_desc.rst'), 'r').read()
# Locally the floe can be invoked by running the terminal command:
# python floes/PlainMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['General MD']]
job.uuid = "266481fc-b257-41e9-b2f9-a92bf028b701"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="solute", title='Solute Input File',
                      description="Solute input file")

sysid = IDSettingCube("System Ids")

md_comp = MDComponentCube("MD Components")
md_comp.set_parameters(multiple_flasks=True)


# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = ParallelSolvationCube("Hydration", title="Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')
solvate.set_parameters(close_solvent=True)

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.2.0')

prod = ParallelMDNptCube("Production", title="Production")
prod.promote_parameter('time', promoted_name='prod_ns', default=2.0,
                       description='Length of MD run in nanoseconds')
prod.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                       description='Temperature (Kelvin)')
prod.promote_parameter('pressure', promoted_name='pressure', default=1.0, description='Pressure (atm)')
prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval', default=0.004,
                       description='Trajectory saving interval in ns')
prod.promote_parameter('hmr', promoted_name="HMR", title='Use Hydrogen Mass Repartitioning', default=False,
                       description='Give hydrogens more mass to speed up the MD')
prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                       description='Select the MD Engine')
prod.set_parameters(reporter_interval=0.004)
prod.set_parameters(suffix='prod')


# Minimization
minComplex = ParallelMDMinimizeCube('minComplex', title='Minimization')
minComplex.set_parameters(restraints="noh (ligand or protein)")
minComplex.set_parameters(restraintWt=5.0)
minComplex.modify_parameter(minComplex.steps, promoted=False, default=0)
minComplex.set_parameters(center=True)
minComplex.set_parameters(save_md_stage=True)
minComplex.set_parameters(hmr=False)
minComplex.promote_parameter("md_engine", promoted_name="md_engine")


# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = ParallelMDNvtCube('warmup', title='System Warm Up')
warmup.set_parameters(time=0.01)
warmup.promote_parameter("temperature", promoted_name="temperature")
warmup.set_parameters(restraints="noh (ligand or protein)")
warmup.set_parameters(restraintWt=2.0)
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
equil1.promote_parameter("hmr", promoted_name="HMR", default=False)
equil1.promote_parameter("temperature", promoted_name="temperature")
equil1.promote_parameter("pressure", promoted_name="pressure")
equil1.set_parameters(restraints="noh (ligand or protein)")
equil1.set_parameters(restraintWt=2.0)
equil1.set_parameters(trajectory_interval=0.0)
equil1.set_parameters(reporter_interval=0.001)
equil1.set_parameters(suffix='equil1')
equil1.promote_parameter("md_engine", promoted_name="md_engine")


# NPT Equilibration stage 2
equil2 = ParallelMDNptCube('equil2', title='Equilibration II')
equil2.set_parameters(time=0.02)
equil2.promote_parameter("temperature", promoted_name="temperature")
equil2.promote_parameter("pressure", promoted_name="pressure")
equil2.promote_parameter("hmr", promoted_name="HMR", default=False)
equil2.set_parameters(restraints="noh (ligand or protein)")
equil2.set_parameters(restraintWt=0.5)
equil2.set_parameters(trajectory_interval=0.0)
equil2.set_parameters(reporter_interval=0.001)
equil2.set_parameters(suffix='equil2')
equil2.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 3
equil3 = ParallelMDNptCube('equil3', title='Equilibration III')
equil3.set_parameters(time=0.03)
equil3.promote_parameter("temperature", promoted_name="temperature")
equil3.promote_parameter("pressure", promoted_name="pressure")
equil3.promote_parameter("hmr", promoted_name="HMR", default=False)
equil3.set_parameters(restraints="ca_protein or (noh ligand)")
equil3.set_parameters(restraintWt=0.1)
equil3.set_parameters(trajectory_interval=0.0)
equil3.set_parameters(reporter_interval=0.001)
equil3.set_parameters(suffix='equil3')
equil3.promote_parameter("md_engine", promoted_name="md_engine")

md_group = ParallelCubeGroup(cubes=[minComplex, warmup, equil1, equil2, equil3, prod])
job.add_group(md_group)

# This cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

rec_check = ParallelRecordSizeCheck("RecordCheck")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(ifs, sysid, md_comp, solvate, coll_open, ff, minComplex,
              warmup, equil1, equil2, equil3, prod,
              coll_close, rec_check, ofs, fail)

ifs.success.connect(sysid.intake)
sysid.success.connect(md_comp.intake)
md_comp.success.connect(solvate.intake)
solvate.success.connect(coll_open.intake)
coll_open.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(prod.intake)
prod.success.connect(coll_close.intake)
coll_close.success.connect(rec_check.intake)
rec_check.success.connect(ofs.intake)

# Fail Connections
sysid.failure.connect(rec_check.fail_in)
solvate.failure.connect(rec_check.fail_in)
md_comp.failure.connect(rec_check.fail_in)
coll_open.failure.connect(rec_check.fail_in)
ff.failure.connect(rec_check.fail_in)
minComplex.failure.connect(rec_check.fail_in)
warmup.failure.connect(rec_check.fail_in)
equil1.failure.connect(rec_check.fail_in)
equil2.failure.connect(rec_check.fail_in)
equil3.failure.connect(rec_check.fail_in)
prod.failure.connect(rec_check.fail_in)
coll_close.failure.connect(rec_check.fail_in)
rec_check.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
