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

from floe.api import WorkFloe

from MDOrion.MDEngines.Gromacs.cubes import (InputGromacs,
                                             GromacsProxyCube,
                                             GromacsRunCube,
                                             WriterRecordCube)

from orionplatform.cubes import DatasetWriterCube

job = WorkFloe('PlainGromacs', title='Plain Gromacs')

job.description = open(path.join(path.dirname(__file__), 'PlainGromacs_desc.rst'), 'r').read()

job.classification = [['General MD']]
job.uuid = "f092b164-7400-403d-8861-b25ff741cab5"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = InputGromacs("Input File", title="Input file")
ifs.promote_parameter('tpr', promoted_name='tpr', default=None)
ifs.promote_parameter("prefix_name", promoted_name="Flask prefix", default="Flask")
ifs.promote_parameter("data_in", promoted_name='in')

proxy = GromacsProxyCube("GromacsProxy", title="Gromacs Proxy Cube")
gmx = GromacsRunCube("GromacsRun", title="Gromacs Run")
gmx.promote_parameter("verbose", promoted_name="verbose", default=False)

ofs = WriterRecordCube("OutputRecords", title="Output Records")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")


job.add_cubes(ifs, proxy, gmx, ofs, fail)

ifs.success.connect(proxy.intake)
proxy.success.connect(gmx.intake)
gmx.success.connect(proxy.intake)
gmx.success.connect(ofs.intake)

# Fail Connections
ifs.failure.connect(fail.intake)
proxy.failure.connect(fail.intake)
gmx.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
