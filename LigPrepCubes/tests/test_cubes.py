import unittest, parmed, base64, pickle
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization, FREDDocking
import OpenMMCubes.utils as utils
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem, oedocking

class ChargeMCMolTester(unittest.TestCase):
    """
    Test the Charging/OMEGA cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = ChargeMCMol('charge')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        ligand = utils.get_data_filename('examples','data/TOL.ism')

        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % ligand)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check outmol has at least 1 conformer
        outmol = self.runner.outputs["success"].get()
        self.assertGreaterEqual(outmol.GetMaxConfIdx(), 1)

        # Loop through atoms and make sure partial charges were set
        for iat, oat in zip(mol.GetAtoms(), outmol.GetAtoms()):
            self.assertNotEqual(iat.GetPartialCharge(), oat.GetPartialCharge)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

class LigandParamTester(unittest.TestCase):
    """
    Test the SMIRNOFF Parameterization cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = LigandParameterization('lig_param')

        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success_gaff(self):
        print('Testing cube:', self.cube.name)
        self.cube.args.molecule_forcefield = 'GAFF'
        self.cube.args.ligand = utils.get_data_filename('examples','data/TOL-chgdmc.oeb.gz')
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(self.cube.args.ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % self.cube.args.ligand)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Get the output molecule
        outmol = self.runner.outputs["success"].get()

        # Check for the attached ParmEd Structure and IDTag
        self.assertTrue(outmol.HasData(oechem.OEGetTag('IDTag')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('Structure')))

        # Check that the Structure is parameterized
        encoded_structure = outmol.GetData(oechem.OEGetTag('Structure'))
        self.assertIsInstance(encoded_structure, str)
        decoded_structure = base64.b64decode(encoded_structure)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        self.assertEqual(len(struct.atoms),mol.NumAtoms())
        struct_string = [x.strip() for x in str(struct).split(';')]
        self.assertIn('parametrized>', struct_string)

    def test_success_smirnoff(self):
        print('Testing cube:', self.cube.name)
        self.cube.args.molecule_forcefield = 'SMIRNOFF'
        self.cube.args.ligand = utils.get_data_filename('examples','data/TOL-chgdmc.oeb.gz')
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(self.cube.args.ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % self.cube.args.ligand)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Get the output molecule
        outmol = self.runner.outputs["success"].get()

        # Check for the attached ParmEd Structure and IDTag
        self.assertTrue(outmol.HasData(oechem.OEGetTag('IDTag')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('Structure')))

        # Check that the Structure is parameterized
        encoded_structure = outmol.GetData(oechem.OEGetTag('Structure'))
        self.assertIsInstance(encoded_structure, str)
        decoded_structure = base64.b64decode(encoded_structure)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        self.assertEqual(len(struct.atoms),mol.NumAtoms())
        struct_string = [x.strip() for x in str(struct).split(';')]
        self.assertIn('parametrized>', struct_string)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class FREDTester(unittest.TestCase):
    """
    Test the FRED docking cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = FREDDocking('fred')
        self.cube.args.receptor = utils.get_data_filename('examples','data/T4-receptor.oeb.gz')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(utils.get_data_filename('examples', 'data/TOL-smnf.oeb.gz'))
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Get the output molecule, check that it has score.
        outmol = self.runner.outputs["success"].get()
        self.assertTrue(oechem.OEHasSDData(outmol, 'Chemgauss4'))

    def test_failure(self):
        pass
    def tearDown(self):
        self.runner.finalize()

if __name__ == "__main__":
        unittest.main()
