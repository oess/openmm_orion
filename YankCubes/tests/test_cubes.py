import unittest, os, parmed
from YankCubes.utils import download_dataset_to_file, get_data_filename
from YankCubes.cubes import YankHydrationCube, YankBindingCube
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem


class YankHydrationCubeTester(unittest.TestCase):
    """
    Test the Yank hydration free energy calculation cube
    """
    def setUp(self):
        self.cube = YankHydrationCube("yank_hydration")
        self.cube.args.timestep = 2 # fs
        self.cube.args.nsteps_per_iteration = 5
        self.cube.args.simulation_time = 0.0001 # ns/replica
        self.cube.args.verbose = True
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        """
        Test that YankHydrationCube can successfully process a single molecule.
        """
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('freesolv_mini.oeb.gz'))

        # Test a single molecule
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()
        # Check that the number of atoms in input and output molecules match.
        self.assertEqual(outmol.NumAtoms(), mol.NumAtoms())
        # Check that a free energy of hydration has been attached
        self.assertTrue(oechem.OEHasSDData(outmol, 'DeltaG_yank_hydration'))
        self.assertTrue(oechem.OEHasSDData(outmol, 'dDeltaG_yank_hydration'))

    def test_failure(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('uranium-hexafluoride.sdf'))

        # Test a single molecule
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 0)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 1)

        outmol = self.runner.outputs["failure"].get()
        # Check that the number of atoms in input and output molecules match.
        self.assertEqual(outmol.NumAtoms(), mol.NumAtoms())
        # Check that an error message has been attached
        self.assertTrue(outmol.HasData(oechem.OEGetTag('error')))

    def test_ports(self):
        pass

    def tearDown(self):
        self.runner.finalize()

class YankBindingCubeTester(unittest.TestCase):
    """
    Test the Yank binding free energy calculation cube
    """
    def setUp(self):
        self.cube = YankBindingCube("yank_binding")
        self.cube.args.receptor = get_data_filename('T4-protein.pdb')
        self.cube.args.timestep = 2 # fs
        self.cube.args.nsteps_per_iteration = 5
        self.cube.args.simulation_time = 0.0001 # ns/replica
        self.cube.args.minimize = False # don't minimize for testing
        self.cube.args.verbose = True
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        """
        Test that YankBindingCube can successfully process a single molecule.
        """
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('p-xylene.mol2'))
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()
        # Check that the number of atoms in input and output molecules match.
        self.assertEqual(outmol.NumAtoms(), mol.NumAtoms())
        # Check that a free energy of hydration has been attached
        self.assertTrue(oechem.OEHasSDData(outmol, 'DeltaG_yank_binding'))
        self.assertTrue(oechem.OEHasSDData(outmol, 'dDeltaG_yank_binding'))

    def test_failure(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('uranium-hexafluoride.sdf'))

        # Test a single molecule
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 0)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 1)

        outmol = self.runner.outputs["failure"].get()
        # Check that the number of atoms in input and output molecules match.
        self.assertEqual(outmol.NumAtoms(), mol.NumAtoms())
        # Check that an error message has been attached
        self.assertTrue(outmol.HasData(oechem.OEGetTag('error')))

    def test_ports(self):
        pass

    def tearDown(self):
        self.runner.finalize()

if __name__ == "__main__":
    unittest.main()
