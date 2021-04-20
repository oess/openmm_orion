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
#
import os

from artemis.wrappers import (WorkFloeWrapper,
                              DatasetWrapper,
                              OutputDatasetWrapper,
                              FileWrapper)

from artemis.test import FloeTestCase

from artemis.decorators import package

import pytest

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

os.chdir(FILE_DIR)


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.local
    @pytest.mark.orion
    def test_omm_PlainMD_protein_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "PlainMD.py"),
            run_timeout=43200,
            queue_timeout=2000
        )
        protein_file = DatasetWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "4JOO_truncNoLig.pdb"
            )
        )

        output = OutputDatasetWrapper(extension=".oedb")
        fail_output = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "solute": protein_file.identifier,
                    "prod_ns": 1,
                    "out": output.identifier,
                    "fail": fail_output.identifier
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # The fail record must be empty
        self.assertEqual(fail_output.count, 0)

        # Check the out record list
        self.assertEqual(output.count, 1)

    @pytest.mark.local
    @pytest.mark.orion
    def test_omm_PlainMD_small_molecule_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "PlainMD.py"),
            run_timeout=43200,
            queue_timeout=2000
        )
        protein_file = DatasetWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "MCL1_lig26.oeb"
            )
        )

        output = OutputDatasetWrapper(extension=".oedb")
        fail_output = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "solute": protein_file.identifier,
                    "prod_ns": 1,
                    "out": output.identifier,
                    "fail": fail_output.identifier
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # The fail record must be empty
        self.assertEqual(fail_output.count, 0)

        # Check the out record list
        self.assertEqual(output.count, 1)

    @pytest.mark.local
    @pytest.mark.orion
    def test_gmx_PlainMD_protein_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "PlainMD.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        protein_file = DatasetWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "4JOO_truncNoLig.pdb"
            )
        )

        output = OutputDatasetWrapper(extension=".oedb")
        fail_output = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "solute": protein_file.identifier,
                    "md_engine": "Gromacs",
                    "prod_ns": 1,
                    "out": output.identifier,
                    "fail": fail_output.identifier
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # The fail record must be empty
        self.assertEqual(fail_output.count, 0)

        # Check the out record list
        self.assertEqual(output.count, 1)

    @pytest.mark.local
    @pytest.mark.orion
    def test_gmx_tpr_floe(self):

        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "PlainGromacs.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        user_tpr_file = FileWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "pL99_ltoluene.tpr"
            )
        )

        fail_output = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "tpr": user_tpr_file.identifier,
                    "fail": fail_output.identifier
                },

                "cube": {
                    "GromacsRun": {
                        "cube_run_time": 0.1
                    }
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # The fail record must be empty
        self.assertEqual(fail_output.count, 0)
