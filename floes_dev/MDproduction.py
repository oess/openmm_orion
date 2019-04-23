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
from MDOrion.MDEngines.cubes import MDNptCube
from cuberecord import DatasetReaderCube, DatasetWriterCube

job = WorkFloe("NPT Production",
               title="NPT Production")

job.description = """
NPT simulation of an OpenMM-ready System

Ex: python floes/openmm_MDnpt.py --system complex.oeb --nanoseconds 0.01

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature in K
pressure (decimal): target final pressure in atm

Outputs:
--------
ofs: Outputs the constant temperature and pressure system
"""

job.classification = [['NPT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

npt = MDNptCube('npt', title="Production")
npt.promote_parameter('time', promoted_name='nanoseconds', default=0.01,
                      description='Length of MD run in nanoseconds')
npt.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                      description='Selected temperature in K')
npt.promote_parameter('pressure', promoted_name='pressure', default=1.0,
                      description='Selected pressure in atm')
npt.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                      description='Select the MD Engine')

# Trajectory and logging info frequency intervals
npt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=0.001,
                      description='Trajectory saving interval in ns')
npt.promote_parameter('reporter_interval', promoted_name='reporter_interval', default=0.001,
                      description='Reporter saving interval in ns')
npt.set_parameters(suffix='prod')

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs, npt, ofs, fail)

ifs.success.connect(npt.intake)
npt.success.connect(ofs.intake)
npt.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
