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

import unittest

import os

from floe.test import CubeTestRunner

import pytest

from MDOrion.MDEngines.cubes import (MDMinimizeCube,
                                     MDNvtCube,
                                     MDNptCube)

import MDOrion

from datarecord import read_records

from openeye import oechem

from MDOrion.Standards import Fields


PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class GmxMinimizationCubeTester(unittest.TestCase):
    """
    Test the Gromacs Minimization cube
    """

    def setUp(self):
        self.cube = MDMinimizeCube('minComplex')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

        os.chdir(FILE_DIR)

    def _test_success(self):
        print('Testing cube:', self.cube.name)

        # Complex file name
        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "pbace_lcat13a.oedb"))

        for record in read_records(ifs):
            pass

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()

        stages = record.get_value(Fields.md_stages)
        self.assertEqual(len(stages), 2)

    @pytest.mark.local
    def test_success(self):
        self.cube.args.md_engine = "Gromacs"
        self.cube.args.steps = 100000
        self._test_success()


class GmxNVTCubeTester(unittest.TestCase):
    """
    Test the Gromacs NVT cube
    """
    def setUp(self):
        self.cube = MDNvtCube('NVT')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

        os.chdir(FILE_DIR)

    def _test_success(self):
        print('Testing cube:', self.cube.name)

        # Complex file name
        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "pP38_lig38a_2n_nvt_5ns.oedb"))

        for record in read_records(ifs):
            pass

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()

        stages = record.get_value(Fields.md_stages)
        self.assertEqual(len(stages), 3)

    @pytest.mark.local
    def test_success(self):
        self.cube.args.md_engine = "Gromacs"
        self.cube.args.time = 0.01  # in nanoseconds
        self.cube.args.nonbondedCutoff = 10.0  # in A
        self.cube.args.temperature = 300.0  # in K
        self.cube.args.restraints = ""
        self.cube.args.save_md_stage = True
        self.cube.args.constraints = "Bonds2H"
        self.cube.args.trajectory_interval = 0.0
        self.cube.args.reporter_interval = 0.0
        self._test_success()


class GmxNPTCubeTester(unittest.TestCase):
    """
    Test the Gromacs NPT cube
    """

    def setUp(self):
        self.cube = MDNptCube('NPT')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

        os.chdir(FILE_DIR)

    def _test_success(self):
        print('Testing cube:', self.cube.name)

        # Complex file name
        ifs = oechem.oeifstream(os.path.join(FILE_DIR, "pP38_lig38a_2n_npt_5ns.oedb"))

        for record in read_records(ifs):
            pass

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()

        stages = record.get_value(Fields.md_stages)
        self.assertEqual(len(stages), 3)

    @pytest.mark.local
    def test_success(self):
        self.cube.args.md_engine = "Gromacs"
        self.cube.args.time = 0.01  # in nanoseconds
        self.cube.args.nonbondedCutoff = 10.0  # in A
        self.cube.args.temperature = 300.0  # in K
        self.cube.args.pressure = 1.0  # in atm
        self.cube.args.restraints = ""
        self.cube.args.save_md_stage = True
        self.cube.args.constraints = "Bonds2H"
        self.cube.args.trajectory_interval = 0.0
        self.cube.args.reporter_interval = 0.0

        self._test_success()


if __name__ == "__main__":
        unittest.main()
