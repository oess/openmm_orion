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


from artemis.wrappers import using_orion

num_proc = 5

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")


os.chdir(FILE_DIR)


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.orion
    def test_EQ_and_NES_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "Equilibrium_and_NES.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "Thrombin3Series_5ligs.oeb"
            )
        )

        protein_file = DatasetWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "Thrombin.pdb"
            )
        )

        map_file = FileWrapper.from_file(
            os.path.join(
                FILE_DIR,
                "thrombin_one_edge.txt"
            )
        )

        output_bound_eq = OutputDatasetWrapper(extension=".oedb")
        output_unbound_eq = OutputDatasetWrapper(extension=".oedb")
        output_nes = OutputDatasetWrapper(extension=".oedb")
        output_fail = OutputDatasetWrapper(extension=".oedb")

        if using_orion():
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "map": map_file.identifier,
                        "out_bound": output_bound_eq.identifier,
                        "out_unbound": output_unbound_eq.identifier,
                        "out": output_nes.identifier,
                        "fail": output_fail.identifier
                    }
                }
            )
        else:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "map": map_file.identifier,
                        "out_bound": output_bound_eq.identifier,
                        "out_unbound": output_unbound_eq.identifier,
                        "out": output_nes.identifier,
                        "fail": output_fail.identifier
                    },

                    "mp": num_proc
                }
            )

        self.assertWorkFloeComplete(workfloe)

        # Check the out record list
        self.assertEqual(output_bound_eq.count, 5)

        # Check the out record list
        self.assertEqual(output_unbound_eq.count, 5)

        # The fail record must be empty
        self.assertEqual(output_fail.count, 0)

        # Check the out record list
        self.assertEqual(output_nes.count, 1)
