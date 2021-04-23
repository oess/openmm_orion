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
                              OutputDatasetWrapper)

from artemis.test import FloeTestCase

from artemis.decorators import package

from openeye.oechem import oeifstream

from datarecord import read_records

import MDOrion

from MDOrion.Standards import Fields

from oeommtools import utils as oeommutils

from openeye import oechem

import pytest

from MDOrion.Standards.mdrecord import MDDataRecord

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DEV_DIR = os.path.join(PACKAGE_DIR, "floes_dev")

os.chdir(FILE_DIR)


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.local
    def test_compex_prep_split_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Complex_prep.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "MCL1_lig26.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "MCL1_protein_ACE_NMA_caps.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "ligands": ligand_file.identifier,
                    "protein": protein_file.identifier,
                    "out": output_file.identifier,
                    "fail": fail_output_file.identifier
                }

            }
        )

        self.assertWorkFloeComplete(workfloe)

        ifs = oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)

        # Each record should have the MD record interface
        for record in records:

            mdrecord = MDDataRecord(record)

            self.assertTrue(record.has_value(Fields.flaskid))
            self.assertTrue(record.has_value(Fields.title))
            self.assertTrue(record.has_value(Fields.ligand))
            self.assertTrue(record.has_value(Fields.protein))
            self.assertTrue(record.has_value(Fields.primary_molecule))
            self.assertTrue(record.has_value(Fields.flask))
            self.assertTrue(record.has_value(Fields.md_stages))
            self.assertTrue(record.has_value(Fields.pmd_structure))

            self.assertEqual(mdrecord.get_flask_id, 0)
            self.assertEqual(mdrecord.get_title, "pMCL1_l26")
            self.assertEqual(record.get_value(Fields.ligand).NumAtoms(), 43)
            self.assertEqual(record.get_value(Fields.protein).NumAtoms(), 2432)

            complx = mdrecord.get_flask
            protein_split, ligand_split, water, excipients = oeommutils.split(complx)
            self.assertEqual(protein_split.NumAtoms(), 2432)
            self.assertEqual(ligand_split.NumAtoms(), 43)
            self.assertEqual(water.NumAtoms(), 20022)
            self.assertEqual(excipients.NumAtoms(), 17)

            stages = mdrecord.get_stages
            self.assertEqual(len(stages), 1)

            stage = stages[0]

            self.assertTrue(stage.has_value(Fields.stage_name))
            self.assertTrue(stage.has_value(Fields.mddata))

            self.assertEqual(stage.get_value(Fields.stage_type), "SETUP")

            top_mol = mdrecord.get_stage_topology()
            self.assertEqual(top_mol.NumAtoms(), complx.NumAtoms())

    @pytest.mark.local
    def test_compex_prep_md_comp_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Complex_prep.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "retigabine.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "kcnq_fix_ARG_fc.oeb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "ligands": ligand_file.identifier,
                    "protein": protein_file.identifier,
                    "out": output_file.identifier,
                    "fail": fail_output_file.identifier
                }

            }
        )

        self.assertWorkFloeComplete(workfloe)

        ifs = oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)

        for record in records:

            md_components = record.get_value(Fields.md_components)

            protein = md_components.get_protein
            ligand = md_components.get_ligand
            counter_ions = md_components.get_counter_ions
            excipients = md_components.get_excipients
            water = md_components.get_water

            self.assertEqual(protein.NumAtoms(), 7353)
            self.assertEqual(ligand.NumAtoms(), 40)
            self.assertEqual(water.NumAtoms(), 58308)
            self.assertEqual(excipients.NumAtoms(), 44)
            self.assertEqual(counter_ions.NumAtoms(), 54)
