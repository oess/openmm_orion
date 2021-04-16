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


import unittest

import MDOrion

from MDOrion.ForceField.cubes import ForceFieldCube

from floe.test import CubeTestRunner

import os

import pytest

from datarecord import read_record

from openeye import oechem

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class ForceFieldPrepTester(unittest.TestCase):
    """
    Test the Complex Preparation cube
    """

    def setUp(self):
        self.cube = ForceFieldCube('ForceFieldPrep')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    @pytest.mark.local
    def test_Gaff2(self):
        print('Testing cube:', self.cube.name)
        # File name
        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'Gaff2'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_Smirnoff99Frosst(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'Smirnoff99Frosst'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_OpenFF1_1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.1.1'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_OpenFF1_2_1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.2.1'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_OpenFF1_3_0(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.0'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_OpenFF1_3_1a1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.1a1'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_protein_non_std_residue(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_protein_force_field_amber_99sbildn_Gaff2(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'Gaff2'
        self.cube.args.protein_forcefield = 'Amber99SBildn'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_99sbildn_OFF_1_2_1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.2.1'
        self.cube.args.protein_forcefield = 'Amber99SBildn'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_99sbildn_OFF_1_3_0(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.0'
        self.cube.args.protein_forcefield = 'Amber99SBildn'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_99sbildn_OFF_1_3_1a1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.1a1'
        self.cube.args.protein_forcefield = 'Amber99SBildn'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    @pytest.mark.local
    def test_protein_force_field_amber_14_Gaff2(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'Gaff2'
        self.cube.args.protein_forcefield = 'Amber14SB'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_14_OFF_1_2_1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.2.1'
        self.cube.args.protein_forcefield = 'Amber14SB'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_14_OFF_1_3_0(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.0'
        self.cube.args.protein_forcefield = 'Amber14SB'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_14_OFF_1_3_1a1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.1a1'
        self.cube.args.protein_forcefield = 'Amber14SB'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_fb15_Gaff2(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'Gaff2'
        self.cube.args.protein_forcefield = 'AmberFB15'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

        @pytest.mark.local
        def test_protein_force_field_amber_fb15_OFF_1_2_1(self):
            print('Testing cube:', self.cube.name)

            ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

            record = read_record(ifs)

            # Selecting ligand and excipient parametrization
            self.cube.args.ligand_forcefield = 'OpenFF_1.2.1'
            self.cube.args.protein_forcefield = 'AmberFB15'

            # Process the molecules
            self.cube.process(record, self.cube.intake.name)

            # Assert that one molecule was emitted on the success port
            self.assertEqual(self.runner.outputs['success'].qsize(), 1)
            # Assert that zero molecules were emitted on the failure port
            self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

            # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_fb15_OFF_1_3_0(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.0'
        self.cube.args.protein_forcefield = 'AmberFB15'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()

    @pytest.mark.local
    def test_protein_force_field_amber_fb15_OFF_1_3_1a1(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "6puq_solvated.oedb"))

        record = read_record(ifs)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'OpenFF_1.3.1a1'
        self.cube.args.protein_forcefield = 'AmberFB15'

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()


if __name__ == "__main__":
    unittest.main()

