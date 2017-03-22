import unittest, os, parmed
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem
import OpenMMCubes.utils as utils

class SetupCubeTester(unittest.TestCase):
    """
    Test the OpenMM complex setup cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = OpenMMComplexSetup("complex_setup")
        self.cube.args.protein = utils.get_data_filename('examples', 'data/epox_hydrolase_apo-protein.pdb')
        self.cube.args.solvent_padding = 1
        self.cube.args.salt_concentration = 10
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(utils.get_data_filename('examples','data/JF6_1-smirff.oeb.gz'))
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check for ParmEd Structure
        outmol = self.runner.outputs["success"].get()
        self.assertTrue(outmol.HasData(oechem.OEGetTag('Structure')))
        # Check structure is for protein:ligand complex, not molecule
        gd = utils.PackageOEMol.unpack(mol)
        self.assertGreater(len(gd['Structure'].atoms), mol.NumAtoms())

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class SimulationCubeTester(unittest.TestCase):
    """
    Test the OpenMM Simulation cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = OpenMMSimulation("md")
        self.runner = CubeTestRunner(self.cube)
        self.cube.args.steps = 3
        self.cube.args.reporter_interval = 1
        self.cube.args.tarxz = False
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def tearDown(self):
        self.runner.finalize()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        mol = oechem.OEMol()
        fname = get_data_filename('examples','data/JF6_1-complex.oeb.gz')
        ifs = oechem.oemolistream(fname)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read complex from %s' % fname)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()
        self.assertTrue(outmol.HasData(oechem.OEGetTag('Structure')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('State')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('Log')))

    def test_failure(self):
        pass

if __name__ == "__main__":
        unittest.main()
