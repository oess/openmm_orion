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


from floe.api import WorkFloe

from cuberecord import (DatasetWriterCube,
                        DatasetReaderCube)

from MDOrion.MDEngines.cubes import (MDMinimizeCube,
                                     MDNvtCube,
                                     MDNptCube)

from MDOrion.System.cubes import SolvationCube

from MDOrion.ForceField.cubes import ForceFieldCube

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting)

job = WorkFloe('Simple MD',
               title='Simple MD')

job.description = """
The Simple MD protocol performs MD simulations given one or more
complete molecular systems as input, each to be treated in its entirety as a solute.
The solute need to have coordinates, all atoms, and correct chemistry.
Each molecular system can have multiple conformers but each conformer will be
run separately as a different solute.
Proteins need to be prepared to an MD standard: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
The input system is solvated and parametrized according to the
selected force fields. A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up and equilibration stages positional harmonic restraints are
applied. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.

Required Input Parameters:
--------------------------
system (file): dataset of prepared ligands posed in the protein active site.

Outputs:
--------
out:  OERecords
"""
# Locally the floe can be invoked by running the terminal command:
# python floes/ShortTrajMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['Complex Setup', 'FrosstMD', 'MD Traj Analysis']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

sysid = IDSettingCube("System Ids")
job.add_cube(sysid)

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = SolvationCube("Hydration", title="System Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')
solvate.set_parameters(close_solvent=True)

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection")
coll_open.set_parameters(open=True)

# Force Field Application
ff = ForceFieldCube("ForceField", title="System Parametrization")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber99SBildn')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='Gaff2')
ff.promote_parameter('other_forcefield', promoted_name='other_ff', default='Gaff2')
ff.set_parameters(lig_res_name='LIG')

prod = MDNptCube("Production", title="Production")
prod.promote_parameter('time', promoted_name='prod_ns', default=2.0,
                       description='Length of MD run in nanoseconds')
prod.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                       description='Temperature (Kelvin)')
prod.promote_parameter('pressure', promoted_name='pressure', default=1.0, description='Pressure (atm)')
prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval', default=0.002,
                       description='Trajectory saving interval in ns')
prod.promote_parameter('hmr', title='Use Hydrogen Mass Repartitioning', default=False,
                       description='Give hydrogens more mass to speed up the MD')
prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                       description='Select the MD Engine')
prod.set_parameters(reporter_interval=0.002)
prod.set_parameters(suffix='prod')


# Minimization
minComplex = MDMinimizeCube('minComplex', title='System Minimization')
minComplex.set_parameters(restraints="noh (ligand or protein)")
minComplex.set_parameters(restraintWt=5.0)
minComplex.set_parameters(steps=0)
minComplex.set_parameters(center=True)
minComplex.set_parameters(save_md_stage=True)
minComplex.promote_parameter("hmr", promoted_name="hmr")
minComplex.promote_parameter("md_engine", promoted_name="md_engine")


# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = MDNvtCube('warmup', title='System Warm Up')
warmup.set_parameters(time=0.01)
warmup.promote_parameter("temperature", promoted_name="temperature")
warmup.set_parameters(restraints="noh (ligand or protein)")
warmup.set_parameters(restraintWt=2.0)
warmup.set_parameters(trajectory_interval=0.0)
warmup.set_parameters(reporter_interval=0.001)
warmup.set_parameters(suffix='warmup')
warmup.promote_parameter("hmr", promoted_name="hmr")
warmup.set_parameters(save_md_stage=True)
warmup.promote_parameter("md_engine", promoted_name="md_engine")


# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1 = MDNptCube('equil1', title='System Equilibration I')
equil1.set_parameters(time=0.01)
equil1.promote_parameter("temperature", promoted_name="temperature")
equil1.promote_parameter("pressure", promoted_name="pressure")
equil1.promote_parameter("hmr", promoted_name="hmr")
equil1.set_parameters(restraints="noh (ligand or protein)")
equil1.set_parameters(restraintWt=2.0)
equil1.set_parameters(trajectory_interval=0.0)
equil1.set_parameters(reporter_interval=0.001)
equil1.set_parameters(suffix='equil1')
equil1.promote_parameter("md_engine", promoted_name="md_engine")


# NPT Equilibration stage 2
equil2 = MDNptCube('equil2', title='System Equilibration II')
equil2.set_parameters(time=0.02)
equil2.promote_parameter("temperature", promoted_name="temperature")
equil2.promote_parameter("pressure", promoted_name="pressure")
equil2.promote_parameter("hmr", promoted_name="hmr")
equil2.set_parameters(restraints="noh (ligand or protein)")
equil2.set_parameters(restraintWt=0.5)
equil2.set_parameters(trajectory_interval=0.0)
equil2.set_parameters(reporter_interval=0.001)
equil2.set_parameters(suffix='equil2')
equil2.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 3
equil3 = MDNptCube('equil3', title='System Equilibration III')
equil3.set_parameters(time=0.03)
equil3.promote_parameter("temperature", promoted_name="temperature")
equil3.promote_parameter("pressure", promoted_name="pressure")
equil3.promote_parameter("hmr", promoted_name="hmr")
equil3.set_parameters(restraints="ca_protein or (noh ligand)")
equil3.set_parameters(restraintWt=0.1)
equil3.set_parameters(trajectory_interval=0.0)
equil3.set_parameters(reporter_interval=0.001)
equil3.set_parameters(suffix='equil3')
equil3.promote_parameter("md_engine", promoted_name="md_engine")

# This cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection")
coll_close.set_parameters(open=False)

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs, solvate, coll_open, ff, minComplex,
              warmup, equil1, equil2, equil3, prod,
              coll_close, ofs, fail)

ifs.success.connect(sysid.intake)
sysid.success.connect(solvate.intake)
solvate.success.connect(coll_open.intake)
coll_open.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(prod.intake)
prod.success.connect(coll_close.intake)
coll_close.success.connect(ofs.intake)
prod.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
