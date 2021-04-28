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

import os

from artemis.wrappers import (WorkFloeWrapper,
                              DatasetWrapper,
                              OutputDatasetWrapper,
                              FileWrapper)

from artemis.test import FloeTestCase

from artemis.decorators import package

import pytest

from openeye.oechem import oeifstream

from datarecord import read_records

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

os.chdir(FILE_DIR)


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.orion
    def test_eq_and_nes_floe(self):

        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "Equilibrium_and_NES.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "Thrombin3Series_5ligs.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "Thrombin.oeb"
            )
        )

        map_file = FileWrapper.get_file(
            os.path.join(
                FILE_DIR,
                "thrombin_one_edge.txt"
            )
        )

        unbound_file = OutputDatasetWrapper(extension=".oedb")
        bound_file = OutputDatasetWrapper(extension=".oedb")

        output_nes_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "ligands": ligand_file.identifier,
                    "protein": protein_file.identifier,
                    "map": map_file.identifier,
                    "out_unbound": unbound_file.identifier,
                    "out_bound": bound_file.identifier,
                    "out": output_nes_file.identifier,
                    "fail": fail_output_file.identifier
                },

            }
        )

        self.assertWorkFloeComplete(workfloe)

        ifs_unbound = oeifstream(unbound_file.path)
        records = []

        for rec in read_records(ifs_unbound):
            records.append(rec)
        ifs_unbound.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 5)

        ifs_bound = oeifstream(bound_file.path)
        records = []

        for rec in read_records(ifs_bound):
            records.append(rec)
        ifs_bound.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 5)

        ifs_nes = oeifstream(output_nes_file.path)
        records = []

        for rec in read_records(ifs_nes):
            records.append(rec)
        ifs_nes.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)
