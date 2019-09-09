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

from MDOrion.System.cubes import CollectionSetting

from MDOrion.TrjAnalysis.cubes import (ParallelTrajToOEMolCube,
                                       ParallelTrajInteractionEnergyCube,
                                       ParallelTrajPBSACube,
                                       ParallelClusterOETrajCube,
                                       ParallelMDTrajAnalysisClusterReport,
                                       MDFloeReportCube,
                                       NMaxWatersLigProt)

job = WorkFloe('MMPBSA', title='MMPBSA')

job.description = """
TESTING
"""
# Locally the floe can be invoked by running the terminal command:
# python floes/ShortTrajMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['Molecular Dynamics']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
isys = DatasetReaderCube("SystemReader", title="System Reader")
isys.promote_parameter("data_in", promoted_name="system",
                       title="System Input File",
                       description="System file name")

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection")
coll_open.set_parameters(open=True)

max_wat = NMaxWatersLigProt("MaxWaters")
max_wat.promote_parameter('explicit_water', promoted_name='explicit_water',
                          default=False,
                          description='Enable MMPBSA calculation with explicit water')

trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube")
IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube")
PBSACube = ParallelTrajPBSACube("TrajPBSACube")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

report = MDFloeReportCube("report", title="Floe Report")

# This cube is necessary dor the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection")
coll_close.set_parameters(open=False)

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(isys,  coll_open, max_wat,
              trajCube, IntECube, PBSACube,
              clusCube, report_gen, report,
              coll_close,  ofs, fail)

isys.success.connect(coll_open.intake)
coll_open.success.connect(max_wat.intake)
max_wat.success.connect(trajCube.intake)
max_wat.failure.connect(fail.intake)
trajCube.success.connect(IntECube.intake)
trajCube.failure.connect(fail.intake)
IntECube.success.connect(PBSACube.intake)
IntECube.failure.connect(fail.intake)
PBSACube.success.connect(clusCube.intake)
PBSACube.failure.connect(fail.intake)
clusCube.success.connect(report_gen.intake)
clusCube.failure.connect(fail.intake)
report_gen.success.connect(report.intake)
report_gen.failure.connect(fail.intake)
report.failure.connect(fail.intake)
report.success.connect(coll_close.intake)
coll_close.success.connect(ofs.intake)
coll_close.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
