# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
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

from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport
from datarecord import (OERecord, Types, OEField)
from openeye import oechem

import tempfile
import MDOrion
import unittest
import unittest.mock as mock
import os
import logging


CLUSTERS_SVG = "pBETA-SECRETA_lHunt13_14_clusters.svg"
HISTRMSD_SVG = "pBETA-SECRETA_lHunt13_14_histRMSD.svg"
LIGINITPOSE_SVG = "pBETA-SECRETA_lHunt13_14_ligInitPose.svg"
TRAJ_SVG = "pBETA-SECRETA_lHunt13_14_traj.svg"

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data", "floe_report")


class TestMDTrajAnalysisFloeReport(unittest.TestCase):
    def read_svg_info_field(self, r, field, filename):
        full_filename = os.path.join(FILE_DIR, filename)
        with open(full_filename, 'r') as f:
            svg_string = f.read()
            r.set_value(field, svg_string)

    def build_traj_record(self):
        oetrajRecord = OERecord()
        self.read_svg_info_field(oetrajRecord, OEField('TrajSVG', Types.String), TRAJ_SVG)

        clusRecord = OERecord()
        self.read_svg_info_field(clusRecord, OEField('HistSVG', Types.String), HISTRMSD_SVG)
        self.read_svg_info_field(clusRecord, OEField('ClusSVG', Types.String), CLUSTERS_SVG)
        self.read_svg_info_field(clusRecord, OEField('rmsdInitPose', Types.String), LIGINITPOSE_SVG)
        clusRecord.set_value(OEField("nFrames", Types.Int), 1)
        clusRecord.set_value(OEField("HDBSCAN_alpha", Types.Float), 1.0)
        clusRecord.set_value(OEField("ClusterMethod", Types.String), 'method')
        clusRecord.set_value(OEField("nClusters", Types.Int), 3)
        clusRecord.set_value(OEField("ClusterCounts", Types.IntVec), [1, 2, 3])

        record_out = OERecord()
        record_out.set_value(OEField('AnalysesDone', Types.StringVec), ["OETraj", "TrajClus"])
        record_out.set_value(OEField('OETraj', Types.Record), oetrajRecord)
        record_out.set_value(OEField('TrajClus', Types.Record), clusRecord)
        record_out.set_value(OEField('Title_OPLMD', Types.String), 'Title')
        m = oechem.OEGraphMol()
        if oechem.OEParseSmiles(m, 'c1ccccc1'):
            record_out.set_value(OEField('Ligand_OPLMD', Types.Chem.Mol), m)
        return record_out

    def test_mdtraj_floe_report(self):
        record = self.build_traj_record()

        self.assertFalse(os.path.exists('md_clus_report.html'), 
          msg="Please delete md_clus_report.html before running")
        floe_report_cube = MDTrajAnalysisClusterReport('floe_report_cube')
        floe_report_cube.begin()

        # Grrr.. not right.. but works
        floe_report_cube.opt['Logger'] = logging
        with mock.patch.object(floe_report_cube.success, 'emit',
                               wraps=floe_report_cube.success.emit) as monkey:
            floe_report_cube.process(record, None)
            monkey.assert_called_with(record)

        floe_report_cube.end()
        self.assertTrue(os.path.exists('md_clus_report.html'))
