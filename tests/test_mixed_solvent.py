# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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
#
import os

from artemis.wrappers import (WorkFloeWrapper,
                              DatasetWrapper,
                              OutputDatasetWrapper,
                              )

from artemis.test import FloeTestCase

from artemis.decorators import package, orion_only

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

os.chdir(FILE_DIR)


@package(PACKAGE_DIR)
class TestMixSolventOrionFloes(FloeTestCase):

    @orion_only
    def test_mix_solvent_md(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "mixed_solvent_md.py"),
            run_timeout=43200,
            queue_timeout=2000
        )
        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "4JOO_truncNoLig.pdb"
            )
        )
        output = OutputDatasetWrapper()

        workfloe.start(
            {
                "promoted": {
                    "protein": protein_file.identifier,
                    "prod_ns": 1,
                    "out": output.identifier,
                    "restraint_to_reference": False,
                },
                "cube": {
                    "Solvation": {"solvents": ["c1ccncc1"]}
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)
        records = list(output.records())
        self.assertEqual(len(records), 1)

    def test_analysis(self):
        from orionclient.types import Dataset
        from orionclient.session import APISession
        dataset = APISession.get_resource(Dataset, 528858)
        from MDOrion.TrjAnalysis.cubes_occupancy import OccupancyCalculator
        cube = OccupancyCalculator()
        cube.begin()
        for rec in dataset.records():
            cube.process(rec, "intake")
        cube.end()