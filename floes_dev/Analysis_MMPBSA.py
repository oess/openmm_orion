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

from MDOrion.MDEngines.cubes import (ParallelMDMinimizeCube,
                                     ParallelMDNvtCube,
                                     ParallelMDNptCube)

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.System.cubes import ParallelSolvationCube

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.ProtPrep.cubes import ProteinSetting

from MDOrion.LigPrep.cubes import (ParallelLigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import (IDSettingCube,
                                  CollectionSetting,
                                  ParallelRecordSizeCheck)

from MDOrion.TrjAnalysis.cubes import (ParallelTrajToOEMolCube,
                                       ParallelTrajInteractionEnergyCube,
                                       ParallelTrajPBSACube,
                                       ParallelClusterOETrajCube,
                                       ParallelMDTrajAnalysisClusterReport,
                                       MDFloeReportCube)

job = WorkFloe('MMPBSA OLD',
               title='MMPBSA OLD')

job.description = """
The Short Trajectory MD (STMD) protocol performs MD simulations given
a prepared protein and a set of posed and prepared ligands as input.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
The protein needs to be prepared to MD standards: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is peformed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up, and equilibration stages, positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.
The production run is then analyzed in terms of interactions between the
ligand and the active site and in terms of ligand RMSD after fitting the trajectory
based on active site C_alphas.

Required Input Parameters:
--------------------------
ligands (file): dataset of prepared ligands posed in the protein active site.
protein (file): dataset of the prepared protein structure.

Outputs:
--------
floe report: html page of the Analysis of each ligand.
out (.oedb file): file of the Analysis results for all ligands.
"""
# Locally the floe can be invoked by running the terminal command:
# python floes/ShortTrajMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['Molecular Dynamics']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
isys = DatasetReaderCube("SystemReader", title="System Reader")
isys.promote_parameter("data_in", promoted_name="system", title="System Input File", description="Ligand file name")

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection")
coll_open.set_parameters(open=True)

trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube")
IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube")
PBSACube = ParallelTrajPBSACube("TrajPBSACube")
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

report = MDFloeReportCube("report", title="Floe Report")

# This cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(isys, coll_open, trajCube, IntECube, PBSACube, clusCube, report_gen,
              report, coll_close, check_rec, ofs, fail)

isys.success.connect(coll_open.intake)
coll_open.success.connect(trajCube.intake)
trajCube.success.connect(IntECube.intake)
trajCube.failure.connect(check_rec.fail_in)
IntECube.success.connect(PBSACube.intake)
IntECube.failure.connect(check_rec.fail_in)
PBSACube.success.connect(clusCube.intake)
PBSACube.failure.connect(check_rec.fail_in)
clusCube.success.connect(report_gen.intake)
clusCube.failure.connect(check_rec.fail_in)
report_gen.success.connect(report.intake)
report_gen.failure.connect(check_rec.fail_in)
report.failure.connect(check_rec.fail_in)
report.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
coll_close.failure.connect(check_rec.fail_in)
check_rec.success.connect(ofs.intake)
check_rec.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
