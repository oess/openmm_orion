import unittest
from LigPrepCubes.cubes import LigandSetChargeCube
import OpenMMCubes.utils as utils
from floe.test import CubeTestRunner
from openeye import oechem

from cuberecord import OEField, OERecord
from cuberecord.constants import DEFAULT_MOL_NAME
from datarecord import Types


class LigChargeTester(unittest.TestCase):
    """
    Test ELF10 charge cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = LigandSetChargeCube('elf10charge')
        self.cube.args.max_conforms = 800
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # File name of a charged ligand
        lig_fname = utils.get_data_filename('examples', 'data/lig_CAT13a_chg.oeb.gz')

        # Read OEMol molecule
        ligand = oechem.OEMol()

        ifs = oechem.oemolistream(lig_fname)
        if not oechem.OEReadMolecule(ifs, ligand):
            raise Exception('Cannot read molecule from %s' % lig_fname)
        ifs.close()

        ligand_copy = ligand.CreateCopy()
        # Set the partial charge to zero
        for at in ligand_copy.GetAtoms():
            at.SetPartialCharge(0.0)

        ligand_record = OERecord()
        field_ligand = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        field_ligand_id = OEField("ID", Types.String)
        ligand_record.set_value(field_ligand, ligand_copy)
        ligand_record.set_value(field_ligand_id, ligand.GetTitle())

        # Process the molecules
        self.cube.process(ligand_record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out ligand
        out_record = self.runner.outputs["success"].get()

        field_out_ligand = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        out_ligand = out_record.get_value(field_out_ligand)

        # Loop through atoms and make sure partial charges were set
        for iat, oat in zip(ligand.GetAtoms(), out_ligand.GetAtoms()):
            self.assertNotEqual(iat.GetPartialCharge(), oat.GetPartialCharge)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


if __name__ == "__main__":
        unittest.main()
