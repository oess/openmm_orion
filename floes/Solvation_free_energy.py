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

from MDOrion.System.cubes import (SolvationCube,
                                  CollectionSetting)

from MDOrion.ForceField.cubes import ForceFieldCube

from MDOrion.LigPrep.cubes import (LigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import IDSettingCube

from MDOrion.Yank.cubes import (YankSolvationFECube,
                                YankProxyCube)

from MDOrion.MDEngines.cubes import (MDMinimizeCube,
                                     MDNvtCube,
                                     MDNptCube)

from MDOrion.TrjAnalysis.cubes import MDFloeReportCube


job = WorkFloe("Solvation Free Energy",
               title="Solvation Free Energy")

job.description = """
The Solvation Free Energy protocol performs Solvation Free Energy Calculations (SFEC) on
a set of input ligands using YANK ( http://getyank.org/latest/ ). The ligands need
to have coordinates, correct chemistry and must be neutral. Each ligand can have multiple
conformers, but each conformer will be prepared and treated as a different ligand.
The ligands are solvated in water (or other solvent or solvent mixture) and parametrized
by the selected force field.
Preceding the SFEC is minimization, warm up, and equilibration in the presence of
positional harmonic restraints. The SFEC is then run by YANK with the selected parameters.
The output floe report contains the Solvation Free Energy values and health checks.

Required Input Parameters:
-----------
ligands: Dataset of the prepared ligands

Outputs:
--------
* out : Dataset of the solvated systems with the calculated solvation free energies
* floe report : An analysis of the results for each ligand
"""

job.classification = [['Solvation Free Energy']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("Ligands", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")
job.add_cube(iligs)


chargelig = LigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Calculate ligand partial charges", default=True)
job.add_cube(chargelig)

ligset = LigandSetting("LigandSetting")
ligset.set_parameters(lig_res_name='LIG')
job.add_cube(ligset)

ligid = IDSettingCube("Ligand Ids")
job.add_cube(ligid)

solvate = SolvationCube("Solvation", title="System Solvation")
solvate.promote_parameter("density", promoted_name="density", title="Solution density in g/ml", default=1.0,
                          description="Solution Density in g/ml")
solvate.promote_parameter("solvents", promoted_name="solvents", title="Solvent components",
                          default='tip3p',
                          description="Comma separated smiles strings of solvent components")
solvate.promote_parameter("molar_fractions", promoted_name="molar_fractions",
                          title="Molar fractions",
                          default='1.0',
                          description="Comma separated strings of solvent molar fractions")
solvate.set_parameters(distance_between_atoms=2.5)
solvate.set_parameters(padding_distance=11.0)

job.add_cube(solvate)

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection")
coll_open.set_parameters(open=True)
job.add_cube(coll_open)

ff = ForceFieldCube("ForceField", title="System Parametrization")
ff.promote_parameter('ligand_forcefield', promoted_name='Ligand ForceField', default='Gaff2')
ff.set_parameters(lig_res_name='LIG')
job.add_cube(ff)

# Add YANK Cube
yank_proxy = YankProxyCube("YankProxy", title="Yank Proxy")
yank_proxy.promote_parameter('iterations', promoted_name='iterations', default=1000,
                             description="Total number of Yank iterations")
job.add_cube(yank_proxy)

# First Yank Cube used to build the UI interface
solvationfe = YankSolvationFECube("SolvationFE", title="Yank Solvation")
solvationfe.promote_parameter('iterations', promoted_name='iterations')
solvationfe.promote_parameter('verbose', promoted_name='verbose', default=False)
# solvationfe.promote_parameter('temperature', promoted_name='temperature', default=300.0,
#                               description='Temperature (Kelvin)')
# solvationfe.promote_parameter('pressure', promoted_name='pressure', default=1.0,
#                               description='Pressure (atm)')
solvationfe.promote_parameter('hmr', promoted_name='hmr', default=False,
                              description='Hydrogen Mass Repartitioning')
solvationfe.set_parameters(lig_res_name='LIG')
job.add_cube(solvationfe)

# Minimization
minimize = MDMinimizeCube("Minimize", title="System Minimization")
minimize.set_parameters(restraints='noh ligand')
minimize.set_parameters(restraintWt=5.0)
minimize.set_parameters(center=True)
minimize.promote_parameter("hmr", promoted_name="hmr")
job.add_cube(minimize)


# NVT Warm-up
warmup = MDNvtCube('warmup', title='System Warm Up')
warmup.set_parameters(time=0.02)
# warmup.promote_parameter("temperature", promoted_name="temperature")
warmup.set_parameters(restraints="noh ligand")
warmup.set_parameters(restraintWt=2.0)
warmup.set_parameters(trajectory_interval=0.0)
warmup.set_parameters(reporter_interval=0.0)
warmup.set_parameters(suffix='warmup')
warmup.promote_parameter("hmr", promoted_name="hmr")
job.add_cube(warmup)


# NPT Equilibration stage
equil = MDNptCube('equil', title='System Equilibration')
equil.set_parameters(time=0.02)
# equil.promote_parameter("temperature", promoted_name="temperature")
# equil.promote_parameter("pressure", promoted_name="pressure")
equil.promote_parameter("hmr", promoted_name="hmr")
equil.set_parameters(restraints="noh ligand")
equil.set_parameters(restraintWt=0.1)
equil.set_parameters(trajectory_interval=0.0)
equil.set_parameters(reporter_interval=0.0)
equil.set_parameters(suffix='equil')
job.add_cube(equil)

report = MDFloeReportCube("report", title="Floe Report")
job.add_cube(report)

# This cube is necessary dor the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection")
coll_close.set_parameters(open=False)
job.add_cube(coll_close)

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")
job.add_cube(ofs)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")
job.add_cube(fail)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(ligset.intake)
ligset.success.connect(ligid.intake)
ligid.success.connect(solvate.intake)
solvate.success.connect(coll_open.intake)
coll_open.success.connect(ff.intake)
ff.success.connect(minimize.intake)
minimize.success.connect(warmup.intake)
warmup.success.connect(equil.intake)
equil.success.connect(yank_proxy.intake)
yank_proxy.success.connect(report.intake)
report.success.connect(coll_close.intake)
report.failure.connect(fail.intake)
coll_close.success.connect(ofs.intake)
yank_proxy.failure.connect(fail.intake)
yank_proxy.cycle_out_port.connect(solvationfe.intake)
solvationfe.success.connect(yank_proxy.cycle_in_port)
solvationfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
