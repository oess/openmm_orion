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

from MDOrion.LigPrep.cubes import ParallelLigandChargeCube

from floe.test import CubeTestRunner

import os

import MDOrion

from MDOrion.Standards import Fields

from openeye import oechem

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class LigChargeTester(unittest.TestCase):
    """
    Test ELF10 charge cube
    """
    def setUp(self):
        self.cube = ParallelLigandChargeCube('elf10charge')
        self.cube.args.max_conforms = 800
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    @pytest.mark.travis
    @pytest.mark.local
    def test_success(self):
        print('Testing cube:', self.cube.name)

        # File name of a charged ligand
        lig_fname = os.path.join(FILE_DIR, "lig_CAT13a_chg.oeb.gz")

        # Read OEMol molecule
        ligand = oechem.OEMol()

        with oechem.oemolistream(lig_fname) as ifs:
            oechem.OEReadMolecule(ifs, ligand)

        ligand_copy = ligand.CreateCopy()
        # Set the partial charge to zero
        for at in ligand_copy.GetAtoms():
            at.SetPartialCharge(0.0)

        ligand_record = OERecord()
        ligand_record.set_value(Fields.primary_molecule, ligand_copy)
        ligand_record.set_value(Fields.wellid, 0)
        ligand_record.set_value(Fields.title, ligand_copy.GetTitle())

        # Process the molecules
        self.cube.process(ligand_record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out ligand
        out_record = self.runner.outputs["success"].get()

        out_ligand = out_record.get_value(Fields.primary_molecule)

        # Loop through atoms and make sure partial charges were set
        for iat, oat in zip(ligand.GetAtoms(), out_ligand.GetAtoms()):
            self.assertNotEqual(iat.GetPartialCharge(), oat.GetPartialCharge)


if __name__ == "__main__":
        unittest.main()
