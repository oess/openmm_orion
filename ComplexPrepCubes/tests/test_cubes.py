import unittest
from ComplexPrepCubes.cubes import ComplexSetPrepCube, ForceFieldSetCube
from openeye import oechem
from OpenMMCubes.cubes import utils as ommutils
from floe.test import CubeTestRunner

from cuberecord import OEField, OERecord
from cuberecord.constants import DEFAULT_MOL_NAME
from datarecord import Types


class ComplexPrepTester(unittest.TestCase):
    """
    Test the Complex Preparation  cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = ComplexSetPrepCube('ComplexPrep')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_protein = ommutils.get_data_filename('examples', 'data/Bace_solvated.oeb.gz')
        fn_ligand = ommutils.get_data_filename('examples', 'data/lig_CAT13a_chg.oeb.gz')

        # Read Protein molecule
        protein = oechem.OEMol()

        with oechem.oemolistream(fn_protein) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        protein_record = OERecord()
        field_protein = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        field_protein_id = OEField("ID", Types.String)
        protein_record.set_value(field_protein, protein)
        protein_record.set_value(field_protein_id, protein.GetTitle())

        # Read Ligand molecule
        ligand = oechem.OEMol()

        with oechem.oemolistream(fn_ligand) as ifs:
            oechem.OEReadMolecule(ifs, ligand)

        ligand_record = OERecord()
        field_ligand = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        field_ligand_id = OEField("ID", Types.String)
        ligand_record.set_value(field_ligand, ligand)
        ligand_record.set_value(field_ligand_id, ligand.GetTitle())

        # Process the molecules
        self.cube.process(protein_record, self.cube.intake.name)
        # Why do I have to manually set these on?
        self.cube.check_protein = True
        self.cube.protein = protein
        self.cube.protein_id = protein.GetTitle()
        self.cube.process(ligand_record, self.cube.protein_port)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        complex_record = self.runner.outputs["success"].get()
        field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        complex = complex_record.get_value(field_complex)

        self.assertEquals(complex.GetMaxAtomIdx(), 52312)

    def tearDown(self):
        self.runner.finalize()

    def test_failure(self):
        pass


class ForceFieldPrepTester(unittest.TestCase):
    """
    Test the Complex Preparation  cube
    Example inputs from `openmm_orion/examples/data`
    """

    def setUp(self):
        self.cube = ForceFieldSetCube('ForceFieldPrep')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_excipient_successGaff2(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_complex = ommutils.get_data_filename('examples',
                                                'data/pbace_lcat13a_solvated_complex.oeb.gz')

        # Read Protein molecule
        complex= oechem.OEMol()

        with oechem.oemolistream(fn_complex) as ifs:
            oechem.OEReadMolecule(ifs, complex)

        complex_record = OERecord()
        field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        field_complex_id = OEField("ID", Types.String)
        complex_record.set_value(field_complex, complex)
        complex_record.set_value(field_complex_id, complex.GetTitle())

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'GAFF2'
        self.cube.args.other_forcefield = 'GAFF2'

        # Process the molecules
        self.cube.process(complex_record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # complex = self.runner.outputs["success"].get()
#
    def test_excipient_successSmirnoff(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_complex = ommutils.get_data_filename('examples',
                                                'data/pbace_lcat13a_solvated_complex.oeb.gz')

        # Read Protein molecule
        complex = oechem.OEMol()

        with oechem.oemolistream(fn_complex) as ifs:
            oechem.OEReadMolecule(ifs, complex)

        complex_record = OERecord()
        field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        field_complex_id = OEField("ID", Types.String)
        complex_record.set_value(field_complex, complex)
        complex_record.set_value(field_complex_id, complex.GetTitle())

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'SMIRNOFF'
        self.cube.args.other_forcefield = 'SMIRNOFF'

        # Process the molecules
        self.cube.process(complex_record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    def test_protein_non_std_residue(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_complex = ommutils.get_data_filename('examples',
                                                'data/pCDK2_l1h1q_solvated_complex.oeb.gz')

        # Read Protein molecule
        complex = oechem.OEMol()

        with oechem.oemolistream(fn_complex) as ifs:
            oechem.OEReadMolecule(ifs, complex)

        complex_record = OERecord()
        field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
        field_complex_id = OEField("ID", Types.String)
        complex_record.set_value(field_complex, complex)
        complex_record.set_value(field_complex_id, complex.GetTitle())

        # Process the molecules
        self.cube.process(complex_record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    def tearDown(self):
        self.runner.finalize()

    def test_failure(self):
        pass


if __name__ == "__main__":
        unittest.main()
