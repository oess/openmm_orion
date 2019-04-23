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

from MDOrion.TrjAnalysis.cubes import TrajToOEMolCube

job = WorkFloe("Testing TrajToOEMol")

job.description = """
Testing Floe
#
Ex. python floes/up.py --in STMD_results.oedb
--out ligands_with_trajOEMol.oedb
#
Parameters:
-----------
in (.oedb file): .oedb file of the MD results
#
Outputs:
--------
ofs (.oedb file): file of the MD results with Traj OEMols
"""
#
ifs = DatasetReaderCube("ifs")
#
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")
#
trajCube = TrajToOEMolCube("TrajToOEMolCube")
#
ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
job.add_cubes(ifs, trajCube, ofs)
#
ifs.success.connect(trajCube.intake)
trajCube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
