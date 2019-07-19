import unittest

import MDOrion

import os

from openeye import oechem

from datarecord import (OERecord,
                        OEField,
                        Types)

from MDOrion.Standards import (Fields,
                               MDFileNames,
                               MDStageTypes)

import pytest

from tempfile import TemporaryDirectory

import tarfile

import pickle


from MDOrion.Standards.mdrecord import MDDataRecord

from datarecord import read_records

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class MDRecordTests(unittest.TestCase):
    """
    Testing MD Record API
    """
    def setUp(self):
        fname = os.path.join(FILE_DIR, "mdrecord.oedb")
        ifs = oechem.oeifstream(fname)
        records = []

        for record in read_records(ifs):
            records.append(record)

        self.assertEqual(len(records), 1)

        self.record = OERecord(records[0])

        self.mdrecord = MDDataRecord(records[0])

        self.cwd = os.getcwd()

        os.chdir(FILE_DIR)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_primary(self):
        mol = self.mdrecord.get_value(Fields.primary_molecule)
        self.assertEqual(mol.NumAtoms(), self.mdrecord.get_primary.NumAtoms())

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_primary(self):
        mol = self.mdrecord.get_value(Fields.primary_molecule)
        self.assertTrue(self.mdrecord.set_primary(mol))

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_well(self):
        mol = self.mdrecord.get_value(Fields.well)
        self.assertEqual(mol.NumAtoms(), self.mdrecord.get_well.NumAtoms())

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_well(self):
        mol = self.mdrecord.get_value(Fields.well)
        self.assertTrue(self.mdrecord.set_well(mol))

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_id(self):
        id = self.mdrecord.get_value(Fields.id)
        self.assertEqual(id, 0)

    @pytest.mark.travis
    @pytest.mark.local
    def test_has_id(self):
        self.assertTrue(self.mdrecord.has_id)

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_id(self):
        self.mdrecord.set_id(5)
        self.assertEqual(self.mdrecord.get_id, 5)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_sys_id(self):
        id = self.mdrecord.get_value(Fields.sysid)
        self.assertEqual(id, 0)

    @pytest.mark.travis
    @pytest.mark.local
    def test_has_sys_id(self):
        self.assertTrue(self.mdrecord.has_sys_id)

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_sys_id(self):
        self.mdrecord.set_sys_id(5)
        self.assertEqual(self.mdrecord.get_sys_id, 5)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_conf_id(self):
        id = self.mdrecord.get_value(Fields.confid)
        self.assertEqual(id, 0)

    @pytest.mark.travis
    @pytest.mark.local
    def test_has_conf_id(self):
        self.assertTrue(self.mdrecord.has_conf_id)

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_conf_id(self):
        self.mdrecord.set_conf_id(5)
        self.assertEqual(self.mdrecord.get_conf_id, 5)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_title(self):
        title = self.mdrecord.get_value(Fields.title)
        self.assertEqual(title, 'pPRT_ltoluene')

    @pytest.mark.travis
    @pytest.mark.local
    def test_has_tile(self):
        self.assertTrue(self.mdrecord.has_title)

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_title(self):
        self.mdrecord.set_title("Pippo")
        self.assertEqual(self.mdrecord.get_title, 'Pippo')

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_last_stage(self):
        last_stage = self.mdrecord.get_last_stage
        self.assertEqual(last_stage.get_value(Fields.stage_name), 'Production')
        self.assertEqual(last_stage.get_value(Fields.stage_type), 'NPT')

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stage_by_name(self):
        last_stage = self.mdrecord.get_stage_by_name()
        self.assertEqual(last_stage.get_value(Fields.stage_name), 'Production')
        self.assertEqual(last_stage.get_value(Fields.stage_type), 'NPT')

        param_stage = self.mdrecord.get_stage_by_name(stg_name='System Parametrization')
        self.assertEqual(param_stage.get_value(Fields.stage_name), 'System Parametrization')
        self.assertEqual(param_stage.get_value(Fields.stage_type), 'SETUP')

        param_stage = self.mdrecord.get_stage_by_name(stg_name='System Minimization')
        self.assertEqual(param_stage.get_value(Fields.stage_name), 'System Minimization')
        self.assertEqual(param_stage.get_value(Fields.stage_type), 'MINIMIZATION')

        with self.assertRaises(ValueError):
            self.mdrecord.get_stage_by_name('Error')

    # @pytest.mark.local
    # def test_delete_stage_by_name(self):
    #     new_record = OERecord(self.record)
    #     new_mdrecord = MDDataRecord(new_record)
    #
    #     new_mdrecord.delete_stage_by_name(stg_name='System Minimization')
    #     self.assertFalse(new_mdrecord.has_stage_name('System Minimization'))
    #     self.assertEqual(len(new_mdrecord.get_stages), 2)

    # TODO THIS IS GOING TO DELETE FILES LOCALLY KEEP DISABLED
    @pytest.mark.travis
    @pytest.mark.local
    def test_has_stage_name(self):
        self.assertTrue(self.mdrecord.has_stage_name('Production'))
        self.assertFalse(self.mdrecord.has_stage_name('Error'))

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stage_by_idx(self):
        with self.assertRaises(ValueError):
            self.mdrecord.get_stage_by_idx(5)
        self.assertEqual(self.mdrecord.get_stage_by_idx(0).get_value(Fields.stage_name), 'System Parametrization')

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stage_state(self):

        last_stage = self.mdrecord.get_last_stage

        mddata_fn = os.path.join(FILE_DIR, last_stage.get_value(Fields.mddata))

        with TemporaryDirectory() as out_directory:

            with tarfile.open(mddata_fn) as tar:
                tar.extractall(path=out_directory)

            state_fn = os.path.join(out_directory, MDFileNames.state)

            with open(state_fn, 'rb') as f:
                md_state = pickle.load(f)

        self.assertEqual(md_state.get_positions(), self.mdrecord.get_stage_state().get_positions())
        self.assertEqual(md_state.get_velocities(), self.mdrecord.get_stage_state().get_velocities())
        self.assertEqual(md_state.get_box_vectors(), self.mdrecord.get_stage_state().get_box_vectors())

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stage_topology(self):

        par_stage = self.mdrecord.get_stage_by_idx(0)

        mddata_fn = os.path.join(FILE_DIR, par_stage.get_value(Fields.mddata))

        with TemporaryDirectory() as out_directory:

            with tarfile.open(mddata_fn) as tar:
                tar.extractall(path=out_directory)

            topology_fn = os.path.join(out_directory, MDFileNames.topology)

            topology_mol = oechem.OEMol()

            with oechem.oemolistream(topology_fn) as ifs:
                oechem.OEReadMolecule(ifs, topology_mol)

        topology = self.mdrecord.get_stage_topology(stg_name='System Parametrization')

        for mol_at, top_at in zip(topology_mol.GetAtoms(), topology.GetAtoms()):
            self.assertEqual(mol_at.GetAtomicNum(), top_at.GetAtomicNum())

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stage_info(self):
        last_stage = self.mdrecord.get_last_stage
        info = last_stage.get_value(Fields.log_data)

        self.assertEqual(info, self.mdrecord.get_stage_info())

        min_stage = self.mdrecord.get_stage_by_name(stg_name='System Minimization')
        info = min_stage.get_value(Fields.log_data)

        self.assertEqual(info, self.mdrecord.get_stage_info(stg_name='System Minimization'))

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stage_trajectory(self):
        self.assertTrue(os.path.isfile(self.mdrecord.get_stage_trajectory()))

    @pytest.mark.travis
    @pytest.mark.local
    def test_add_new_stage(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)

        topology = self.mdrecord.get_stage_topology()
        md_state = self.mdrecord.get_stage_state()

        self.assertTrue(new_mdrecord.add_new_stage("Testing",
                                                   MDStageTypes.FEC,
                                                   topology,
                                                   md_state,
                                                   "test.tar.gz",
                                                   log='TestingLogs'))

        self.assertEqual(len(new_mdrecord.get_value(Fields.md_stages)), 4)
        new_last_stage = new_mdrecord.get_stage_by_name(stg_name='Testing')

        self.assertEqual(new_last_stage.get_value(Fields.stage_name), 'Testing')
        self.assertEqual(new_last_stage.get_value(Fields.stage_type), MDStageTypes.FEC)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stages(self):
        stages = self.mdrecord.get_stages
        self.assertEqual(len(stages), 3)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_stages_names(self):
        ls_names = ['System Parametrization', 'System Minimization', 'Production']
        stg_names = self.mdrecord.get_stages_names
        self.assertEqual(stg_names, ls_names)

    @pytest.mark.travis
    @pytest.mark.local
    def test_has_stages(self):
        self.assertTrue(self.mdrecord.has_stages)

    @pytest.mark.travis
    @pytest.mark.local
    def test_get_parmed(self):
        pmd = self.mdrecord.get_parmed(sync_stage_name='last')
        self.assertEqual(len(pmd.atoms), 30439)
        self.assertEqual((len(pmd.residues)), 9446)
        self.assertEqual((len(pmd.bonds)), 21178)
        self.assertEqual((len(pmd.angles)), 14069)
        self.assertEqual((len(pmd.dihedrals)), 8028)

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_parmed(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)

        pmd = new_mdrecord.get_parmed()

        new_mdrecord.delete_field(Fields.pmd_structure)
        self.assertFalse(new_mdrecord.has_parmed)

        new_mdrecord.set_parmed(pmd, sync_stage_name='last')
        self.assertTrue(new_mdrecord.has_parmed)

    @pytest.mark.travis
    @pytest.mark.local
    def test_has_parmed(self):
        self.assertTrue(self.mdrecord.has_parmed)

    @pytest.mark.travis
    @pytest.mark.local
    def test_delete_parmed(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)

        new_mdrecord.delete_parmed
        self.assertFalse(new_mdrecord.has_parmed)

    @pytest.mark.travis
    @pytest.mark.local
    def test_protein_traj(self):

        oetraj_record = self.record.get_value(OEField('OETraj', Types.Record))

        prot_mol = oetraj_record.get_value(Fields.protein_traj_confs)

        mdrecord = MDDataRecord(oetraj_record)

        mdprot = mdrecord.get_protein_traj

        self.assertEqual(prot_mol.NumAtoms(), mdprot.NumAtoms())

    @pytest.mark.travis
    @pytest.mark.local
    def test_set_protein_traj(self):
        oetraj_record = self.record.get_value(OEField('OETraj', Types.Record))

        prot_mol = oetraj_record.get_value(Fields.protein_traj_confs)

        mdrecord = MDDataRecord(oetraj_record)

        oetraj_record.delete_field(Fields.protein_traj_confs)

        with self.assertRaises(ValueError):
            mdrecord.get_protein_traj

        self.assertTrue(mdrecord.set_protein_traj(prot_mol))
