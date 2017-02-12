import unittest, parmed, base64, pickle
from LigPrepCubes.ports import (CustomMoleculeInputPort, CustomMoleculeOutputPort)
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem, oedocking

class OMEGATester(unittest.TestCase):
    """
    Test the OMEGA cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = OEOmegaConfGen('omega')
        self.cube.args.maxConfs = 1
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        self.cube.args.ligand = get_data_filename('test_smiles.ism')

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
        self.assertEqual(outmol.GetMaxConfIdx(), 1)

    def test_failure(self):
        pass
    def tearDown(self):
        self.runner.finalize()

class FREDTester(unittest.TestCase):
    """
    Test the OMEGA cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = FREDDocking('fred')
        self.cube.args.receptor = get_data_filename('test-receptor.oeb.gz')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('mcmol.oeb'))
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Get the output molecule
        outmol = self.runner.outputs["success"].get()
        self.assertTrue(oechem.OEHasSDData(outmol))

    def test_failure(self):
        pass
    def tearDown(self):
        self.runner.finalize()

class SMIRFFTester(unittest.TestCase):
    """
    Test the SMIRFF Parameterization cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = SMIRFFParameterization('smirff')
        self.cube.args.molecule_forcefield = 'smirff99Frosst.ffxml'
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        self.cube.args.ligand = get_data_filename('toluene.pdb')
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

        # Check for the OpenMM System
        serialized_system = outmol.GetData(oechem.OEGetTag('system'))
        # System should be encoded to str on output
        self.assertIsInstance(serialized_system, str)
        # Check it can regenerate System
        system = openmm.XmlSerializer.deserialize(serialized_system)
        self.assertIsInstance(system, openmm.System)
        self.assertEqual(system.getNumParticles(),mol.NumAtoms())

        # Check for the ParmEd Structure
        encoded_structure = outmol.GetData(oechem.OEGetTag('structure'))
        self.assertIsInstance(encoded_structure, str)
        decoded_structure = base64.b64decode(encoded_structure)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        self.assertEqual(len(struct.atoms),mol.NumAtoms())
        print(struct)
        struct_string = [x.strip() for x in str(struct).split(';')]
        self.assertIn('parametrized>', struct_string)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

if __name__ == "__main__":
        unittest.main()
