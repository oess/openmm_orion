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

from openeye import oechem

from datarecord import read_records

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

        output_bound_eq_file = OutputDatasetWrapper(extension=".oedb")
        output_unbound_eq_file = OutputDatasetWrapper(extension=".oedb")
        output_nes_file = OutputDatasetWrapper(extension=".oedb")
        output_fail_file = OutputDatasetWrapper(extension=".oedb")

        if using_orion():
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "map": map_file.identifier,
                        "out_bound": output_bound_eq_file.identifier,
                        "out_unbound": output_unbound_eq_file.identifier,
                        "out": output_nes_file.identifier,
                        "fail": output_fail_file.identifier
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
                        "out_bound": output_bound_eq_file.identifier,
                        "out_unbound": output_unbound_eq_file.identifier,
                        "out": output_nes_file.identifier,
                        "fail": output_fail_file.identifier
                    },

                    "mp": num_proc
                }
            )

        self.assertWorkFloeComplete(workfloe)

        # Check output
        ifs = oechem.oeifstream(output_bound_eq_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 5)

        # Check output
        ifs = oechem.oeifstream(output_unbound_eq_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 5)

        fail_ifs = oechem.oeifstream(output_fail_file.path)
        records_fail = []

        for rec_fail in read_records(fail_ifs):
            records_fail.append(rec_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        # Check output
        ifs = oechem.oeifstream(output_nes_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)
