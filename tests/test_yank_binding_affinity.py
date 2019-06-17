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
                              FileWrapper,
                              OutputDatasetWrapper)

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
FLOES_DEV_DIR = os.path.join(PACKAGE_DIR, "floes_dev")

os.chdir(FILE_DIR)


@package(PACKAGE_DIR)
class TestYankBindingFloes(FloeTestCase):

    @pytest.mark.local
    def test_yank_binding_repex_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Binding_free_energy_repex.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "toluene.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "lysozyme.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        if using_orion:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }

                }
            )
        else:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "mp": num_proc,

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }

                }
            )

        self.assertWorkFloeComplete(workfloe)

        fail_ifs = oechem.oeifstream()
        records_fail = []

        for rec_fail in read_records(fail_ifs):
            records_fail.append(rec_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        # Check output
        ifs = oechem.oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)

    @pytest.mark.local
    def test_yank_binding_repex_multi_ligs_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Binding_free_energy_repex.py"),
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
                "Thrombin.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        if using_orion:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }
                }
            )
        else:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "mp": num_proc,

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }
                }
            )

        self.assertWorkFloeComplete(workfloe)

        fail_ifs = oechem.oeifstream()
        records_fail = []

        for rec_fail in read_records(fail_ifs):
            records_fail.append(rec_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        # Check output
        ifs = oechem.oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 5)

    @pytest.mark.local
    def test_yank_binding_sams_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Binding_free_energy_sams.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "toluene.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "lysozyme.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        if using_orion:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }
                }
            )
        else:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "mp": num_proc,

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }
                }
            )

        self.assertWorkFloeComplete(workfloe)

        fail_ifs = oechem.oeifstream()
        records_fail = []

        for rec_fail in read_records(fail_ifs):
            records_fail.append(rec_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        # Check output
        ifs = oechem.oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)

    @pytest.mark.multi
    @pytest.mark.local
    def test_yank_binding_sams_multi_ligs_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Binding_free_energy_sams.py"),
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
                "Thrombin.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        if using_orion:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }
                }
            )
        else:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "mp": num_proc,

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }
                }
            )

        self.assertWorkFloeComplete(workfloe)

        fail_ifs = oechem.oeifstream()
        records_fail = []

        for rec_fail in read_records(fail_ifs):
            records_fail.append(rec_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        # Check output
        ifs = oechem.oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 5)

    @pytest.mark.local
    def test_yank_binding_yaml_user_file(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "Binding_free_energy_detailed.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "toluene.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "lysozyme.pdb"
            )
        )

        user_yaml_file = FileWrapper.get_file(
            os.path.join(
                FILE_DIR,
                "boresch.yaml"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        if using_orion:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "yaml": user_yaml_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },
                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }

                }
            )
        else:
            workfloe.start(
                {
                    "promoted": {
                        "ligands": ligand_file.identifier,
                        "protein": protein_file.identifier,
                        "yaml": user_yaml_file.identifier,
                        "iterations": 13,
                        "out": output_file.identifier,
                        "fail": fail_output_file.identifier
                    },

                    "mp": num_proc,

                    "cube": {
                        "YankABFE": {
                            "checkpoint_interval": 10
                        }
                    }

                }
            )

        self.assertWorkFloeComplete(workfloe)

        fail_ifs = oechem.oeifstream()
        records_fail = []

        for rec_fail in read_records(fail_ifs):
            records_fail.append(rec_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        # Check output
        ifs = oechem.oeifstream(output_file.path)
        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)
