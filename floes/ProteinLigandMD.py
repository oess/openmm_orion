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

from floes.SubfloeFunctions import setup_MD_startup

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

job = WorkFloe('Solvate and Run Protein-Ligand MD', title='Solvate and Run Protein-Ligand MD')

job.description = open(path.join(path.dirname(__file__), 'ProteinLigandMD_desc.rst'), 'r').read()

job.classification = [['Specialized MD']]
job.uuid = "ae561d76-a2b6-4d89-b621-b979f1930b40"
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

# This Cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)
coll_open.set_parameters(write_new_collection='MD_OPLMD')

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.0')

# Protein Setting
mdcomp = MDComponentCube("MD Components", title="MD Components")
mdcomp.promote_parameter("flask_title", promoted_name="flask_title", default="")

# This Cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(iligs, ligset, iprot, mdcomp, chargelig, complx,
              solvate, coll_open, ff,
              coll_close, check_rec, ofs, fail)

MD_startup_options = {}
MD_startup_options['Prod_Default_Time_ns'] = 2.0
MD_startup_options['Prod_Default_Traj_Intvl_ns'] = 0.004
setup_MD_startup(job, ff, coll_close, check_rec, MD_startup_options)

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
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
