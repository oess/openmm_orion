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

from datarecord import OERecord

import unittest

from MDOrion.ProtPrep.cubes import MDSetting

from floe.test import CubeTestRunner

import os

import MDOrion

from datarecord import read_records

from MDOrion.Standards import Fields

from openeye import oechem

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class MDSettingTester(unittest.TestCase):
    """
    MD Components testing
    """
    def setUp(self):
        self.cube = MDSetting("MDSetting")
        self.cube.args.multiple_flask = False
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    @pytest.mark.travis
    @pytest.mark.local
    def test_success_from_du(self):
        print('Testing cube:', self.cube.name)

        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "1H1Q_du.oedb"))
        for record in read_records(ifs):
            pass
        ifs.close()

        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()
        md_components = record.get_value(Fields.md_components)

        self.assertEqual(md_components.num_atoms, 9746)
        self.assertTrue(md_components.has_protein)
        self.assertTrue(md_components.has_ligand)
        self.assertTrue(md_components.has_solvent)
        self.assertFalse(md_components.has_cofactors)

    @pytest.mark.travis
    @pytest.mark.local
    def test_success_from_molecule_spruce(self):
        print('Testing cube:', self.cube.name)

        protein_fn = os.path.join(FILE_DIR, "4YFF_prot.oeb")
        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        record = OERecord()
        record.set_value(Fields.primary_molecule, protein)

        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()
        md_components = record.get_value(Fields.md_components)

        self.assertEqual(md_components.num_atoms, 8442)
        self.assertTrue(md_components.has_protein)
        self.assertTrue(md_components.has_ligand)
        self.assertTrue(md_components.has_solvent)
        self.assertFalse(md_components.has_cofactors)

    @pytest.mark.travis
    @pytest.mark.local
    def test_success_from_molecule_splitter(self):
        print('Testing cube:', self.cube.name)

        protein_fn = os.path.join(FILE_DIR, "lysozyme.pdb")
        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        record = OERecord()
        record.set_value(Fields.primary_molecule, protein)

        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()
        md_components = record.get_value(Fields.md_components)

        self.assertEqual(md_components.num_atoms, 2635)
        self.assertTrue(md_components.has_protein)
        self.assertFalse(md_components.has_ligand)
        self.assertFalse(md_components.has_solvent)
        self.assertFalse(md_components.has_cofactors)


if __name__ == "__main__":
        unittest.main()


