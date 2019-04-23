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
from MDOrion.MDEngines.cubes import MDMinimizeCube
from cuberecord import DatasetReaderCube, DatasetWriterCube

job = WorkFloe("Minimize",
               title="Minimize")

job.description = """
Minimize an OpenMM-ready solvated complex

Ex: python floes/openmm_prepMDminimize.py --system complex.oeb --ofs-data_out min.oeb --steps 1000`

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
steps (int): Number of MD steps to minimize the system. If 0 until convergence will be reached

Outputs:
--------
ofs: Outputs the minimized system
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

min = MDMinimizeCube('Minimize', title="System Minimization")
min.promote_parameter('steps', promoted_name='steps', default=0)
min.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                      description='Select the MD Engine')
min.set_parameters(save_md_stage=True)

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs, min, ofs, fail)
ifs.success.connect(min.intake)
min.success.connect(ofs.intake)
min.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
