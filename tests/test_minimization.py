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

import pytest

from openeye.oechem import oeifstream

from datarecord import read_records

import MDOrion

from simtk import (unit,
                   openmm)

from simtk.openmm import app

from openeye import oechem

from MDOrion.Standards.mdrecord import MDDataRecord

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DEV_DIR = os.path.join(PACKAGE_DIR, "floes_dev")


os.chdir(FILE_DIR)


# Supporting functions
def calculate_eng(mdstate, parmed_structure):

        topology = parmed_structure.topology
        positions = mdstate.get_positions()
        box = mdstate.get_box_vectors()

        # OpenMM system
        system = parmed_structure.createSystem(nonbondedMethod=app.PME,
                                               nonbondedCutoff=10 * unit.angstroms,
                                               constraints=app.HBonds)
        # OpenMM Integrator
        integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin,
                                               1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
        # Set Simulation
        simulation = app.Simulation(topology, system, integrator)

        # Set Positions
        simulation.context.setPositions(positions)

        # Set Box dimensions
        simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

        # Collect the OpenMM state energy info
        state = simulation.context.getState(getEnergy=True)

        # Potential Energy
        peng = state.getPotentialEnergy()

        return peng


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.local
    def test_omm_minimization_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "MDminimize.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        system = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "pbace_lcat13a.oedb"
            )
        )

        # Read input record
        ifs = oeifstream(system.dataset_path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)

        # Check the out record list
        self.assertEqual(count, 1)

        # Calculate the initial potential energy
        for record in records:

            mdrecord = MDDataRecord(record)

            stages = mdrecord.get_stages
            self.assertEqual(len(stages), 1)

            mdstate = mdrecord.get_stage_state()
            parmed_structure = mdrecord.get_parmed(sync_stage_name='last')

            # Calculate the initial potential energy
            eng_i = calculate_eng(mdstate, parmed_structure)

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "system": system.identifier,
                    "out": output_file.identifier,
                    "fail": fail_output_file.identifier
                },

                "cube": {
                    "Minimize": {
                        "save_md_stage": True
                    }
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # Read output record
        ifs = oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # The records list must have just one record
        self.assertEqual(count, 1)

        # Calculate the final potential energy
        for record in records:

            mdrecord = MDDataRecord(record)

            stages = mdrecord.get_stages
            self.assertEqual(len(stages), 2)

            mdstate = mdrecord.get_stage_state()
            parmed_structure = mdrecord.get_parmed(sync_stage_name='last')

            # Calculate the final potential energy
            eng_f = calculate_eng(mdstate, parmed_structure)

        self.assertLess(eng_f.in_units_of(unit.kilojoule_per_mole)/unit.kilojoule_per_mole,
                        eng_i.in_units_of(unit.kilojoule_per_mole)/unit.kilojoule_per_mole)

    @pytest.mark.local
    def test_gmx_minimization_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "MDminimize.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        system = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "pbace_lcat13a.oedb"
            )
        )

        # Read input record
        ifs = oeifstream(system.dataset_path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)

        # Check the out record list
        self.assertEqual(count, 1)

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "system": system.identifier,
                    "md_engine": "Gromacs",
                    "out": output_file.identifier,
                    "fail": fail_output_file.identifier
                },

                "cube": {
                    "Minimize": {
                        "save_md_stage": True
                    }
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # Read output record
        ifs = oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # The records list must have just one record
        self.assertEqual(count, 1)
