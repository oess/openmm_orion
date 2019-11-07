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

from MDOrion.MDEngines.Gromacs.cubes import (InputGromacs,
                                             GromacsProxyCube,
                                             GromacsRunCube,
                                             WriterRecordCube)


job = WorkFloe('PlainGromacs', title='Plain Gromacs')

job.description = """
This Floe has been design to run Gromacs by using as
input Gromacs .tpr files. The floe will run Gromacs
for a max number of hours (10hrs default) outputting the produced 
trajectory file and a recovery dataset. The Gromacs cube 
runs in a cycle till the number of md steps specified in 
the .tpr file are consumed. If the recovery dataset is 
provided as input, Gromacs will recover and run from the last
check point saved in the recovery dataset. If both .tpr 
and recovery dataset files are provided the recovery dataset 
will overwrite the .tpr file.

Required Input Parameters:
--------------------------
tpr (file): Gromacs Tpr file
dataset: The recovery Gromacs dataset

Outputs:
--------

OEDataset Recovery file
OEFile Gromacs Trajectory files
OEFile Restart file (useless for the user)
"""

job.classification = [['Molecular Dynamics']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = InputGromacs("Input File", title="Input file")
ifs.promote_parameter('tpr', promoted_name='tpr', default=None)
ifs.promote_parameter("prefix_name", promoted_name="prefix", default="PROT")
ifs.promote_parameter("data_in", promoted_name='in')

proxy = GromacsProxyCube("GromacsProxy", title="Gromacs Proxy Cube")
gmx = GromacsRunCube("GromacsRun", title="Gromacs Run")
gmx.promote_parameter("verbose", promoted_name="verbose", default=False)

ofs = WriterRecordCube("OutputRecords", title="Output Records")

job.add_cubes(ifs, proxy, gmx, ofs)

ifs.success.connect(proxy.intake)
proxy.success.connect(gmx.intake)
gmx.success.connect(proxy.intake)
gmx.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
