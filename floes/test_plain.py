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

from MDOrion.System.cubes import ParallelSolvationCube

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.System.cubes import MDComponentCube


from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelRecordSizeCheck)

job = WorkFloe('TESTING',
               title='TESTING')

job.classification = [['Specialized MD']]
job.uuid = "5e0de570-dc8f-4b26-a8f2-c189ba22a56d"
job.tags = [tag for lists in job.classification for tag in lists]


# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
flask = DatasetReaderCube("FaskReader", title="Flask Reader")
flask.promote_parameter("data_in", promoted_name="flask", title='Flask Input Dataset',
                        description="Flask Dataset")

idsetting = IDSettingCube("IDSetting")

# Protein Setting
comp = MDComponentCube("ComponentSetting", title="Component Setting")
comp.promote_parameter("flask_title", promoted_name="protein_title", default="")
comp.set_parameters(multiple_flasks=True)

solvate = ParallelSolvationCube("Solvation", title="Solvation")
solvate.set_parameters(density=1.03)
solvate.set_parameters(salt_concentration=50.0)
solvate.set_parameters(neutralize_solute=True)
solvate.modify_parameter(solvate.close_solvent, promoted=False, default=False)

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.0.0')

minflask = ParallelMDMinimizeCube("Minimization")
minflask.modify_parameter(minflask.restraints, promoted=False, default="noh (ligand or protein)")
minflask.modify_parameter(minflask.restraintWt, promoted=False, default=5.0)
minflask.set_parameters(center=True)
minflask.set_parameters(save_md_stage=True)
minflask.promote_parameter("md_engine", promoted_name="md_engine")


ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(flask, idsetting, comp, solvate, ff, minflask, ofs, fail)

flask.success.connect(idsetting.intake)
idsetting.success.connect(comp.intake)
comp.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(minflask.intake)
minflask.success.connect(ofs.intake)

# Fail Connections
idsetting.failure.connect(fail.intake)
comp.failure.connect(fail.intake)
solvate.failure.connect(fail.intake)
ff.failure.connect(fail.intake)
minflask.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
