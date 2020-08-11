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

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.ForceField.cubes import ParallelEnergyCube

from MDOrion.System.cubes import MDComponentCube


job = WorkFloe('Test Amber FF',
               title='Test Amber FF')

job.classification = [['General MD']]
job.uuid = "054ef9ab-f93a-4a8f-8679-cbf90bfc9def"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                      description="Protein input file")

md_comp = MDComponentCube("MD Components")
md_comp.set_parameters(multiple_flasks=True)

eng = ParallelEnergyCube("EnergyDecomposition")
eng.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(ifs, md_comp, eng, ofs, fail)

ifs.success.connect(md_comp.intake)
md_comp.success.connect(eng.intake)
eng.success.connect(ofs.intake)

md_comp.failure.connect(fail.intake)
eng.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
