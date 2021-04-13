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

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import ExtractMDDataCube

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from os import path

job = WorkFloe("Extract Short Trajectory MD Results for Download",
               title="Extract Short Trajectory MD Results for Download")


job.description = open(path.join(path.dirname(__file__), 'MDData_desc.rst'), 'r').read()

job.classification = [['MD Data']]
job.uuid = "6665ca20-6014-4f3b-8d02-4b5d15b75ee3"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system",
                      title='STMDA Input File',
                      description="The Dataset produced by the Short Trajectory MD with Analysis floe", order=0)

data = ExtractMDDataCube("MDData", title="Extract MD Data")

data.promote_parameter('out_file_name', promoted_name='out_file_name',
                       description="Output File name",
                       default="md_data.tar.gz", oeder=1)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", description="Fail Data Set", order=2)

job.add_cubes(ifs, data, fail)
ifs.success.connect(data.intake)
data.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
