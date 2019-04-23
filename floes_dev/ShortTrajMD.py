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

from cuberecord import ( DatasetWriterCube,
                         DatasetReaderCube)

from MDOrion.MDEngines.cubes import (MDMinimizeCube,
                                     MDNvtCube,
                                     MDNptCube)

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.System.cubes import SolvationCube

from MDOrion.ForceField.cubes import ForceFieldCube

from MDOrion.ProtPrep.cubes import ProteinSetting

from MDOrion.LigPrep.cubes import (LigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import IDSettingCube

from MDOrion.System.cubes import CollectionSetting


job = WorkFloe('Short Trajectory MD',
               title='Short Trajectory MD')

job.description = """
NOTE: this is an Alpha Test version.
We are actively working on improving the MD sampling.

The Short Trajectory MD (STMD) protocol performs MD simulations given a set of
prepared ligands and a prepared protein as input.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
The protein needs to be prepared to an MD standard: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
The STMD floe requires as inputs the protein and a set of ligands correctly posed
in the protein binding site. A complex is formed with each ligand and conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is peformed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up and equilibration stages positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.

Required Input Parameters:
--------------------------
ligands (file): dataset of prepared ligands posed in the protein active site.
protein (file): dataset of the prepared protein structure.

Outputs:
--------
out:  OERecords (one per ligand) of MD and Analysis results.
floe report: html page of the Analysis of each ligand.
"""
# Locally the floe can be invoked by running the terminal command:
# python floes/ShortTrajMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['Complex Setup', 'FrosstMD']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)

ligset = LigandSetting("LigandSetting", title="Ligand Setting")
ligset.set_parameters(lig_res_name='LIG')

ligid = IDSettingCube("Ligand Ids")
job.add_cube(ligid)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = DatasetReaderCube("ProteinReader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                        description="Protein file name")

protset = ProteinSetting("ProteinSetting", title="Protein Setting")
protset.promote_parameter("protein_prefix", promoted_name="protein_prefix", default="PRT")

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex", title="Complex Preparation")
complx.set_parameters(lig_res_name='LIG')

# The solvation cube is used to solvate the system and define the ionic strength of the solution
# solvate = HydrationCube("Hydration")

solvate = SolvationCube("Hydration", title="System Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')
solvate.set_parameters(close_solvent=True)

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
prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                       description='Select the MD Engine')
prod.set_parameters(reporter_interval=0.002)
prod.set_parameters(suffix='prod')


# Minimization
minComplex = MDMinimizeCube('minComplex', title='System Minimization')
minComplex.promote_parameter("hmr", promoted_name="hmr")
minComplex.set_parameters(restraints="noh (ligand or protein)")
minComplex.set_parameters(restraintWt=5.0)
minComplex.set_parameters(steps=0)
minComplex.set_parameters(center=True)
minComplex.set_parameters(save_md_stage=True)
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

coll_close = CollectionSetting("CloseCollection")
coll_close.set_parameters(open=False)

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(iligs, ligset, iprot, protset, chargelig, complx, solvate, coll_open, ff,
              minComplex, warmup, equil1, equil2, equil3, prod, coll_close, ofs, fail)


iligs.success.connect(chargelig.intake)
chargelig.success.connect(ligset.intake)
ligset.success.connect(ligid.intake)
ligid.success.connect(complx.intake)
iprot.success.connect(protset.intake)
protset.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)
solvate.success.connect(coll_open.intake)
coll_open.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
ff.failure.connect(fail.intake)
minComplex.success.connect(warmup.intake)
minComplex.failure.connect(fail.intake)
warmup.success.connect(equil1.intake)
warmup.failure.connect(fail.intake)
equil1.success.connect(equil2.intake)
equil1.failure.connect(fail.intake)
equil2.success.connect(equil3.intake)
equil2.failure.connect(fail.intake)
equil3.success.connect(prod.intake)
equil3.failure.connect(fail.intake)
prod.success.connect(coll_close.intake)
prod.failure.connect(fail.intake)
coll_close.success.connect(ofs.intake)
coll_close.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
