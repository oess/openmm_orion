import unittest, os, parmed
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem


class SetupCubeTester(unittest.TestCase):
    """
    Test the OpenMM complex setup cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = OpenMMComplexSetup("complex_setup")
        self.cube.args.protein = get_data_filename('T4-protein.pdb')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('9PC1X-smirff.oeb.gz'))
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
        self.assertGreater(outmol.NumAtoms(), mol.NumAtoms())
        self.assertTrue(outmol.HasData(oechem.OEGetTag('structure')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('system')))

    def test_failure(self):
        pass

    def test_ports(self):
        mol = oechem.OEMol()
        system = openmm.System()
        structure = parmed.structure.Structure()
        sys_out = OpenMMSystemOutput("sys_out")
        struct_out = ParmEdStructureOutput('struct_out')
        with oechem.oemolostream('empty-mol.oeb.gz') as ofs:
            mol.SetData(oechem.OEGetTag('system'), sys_out.encode(system))
            mol.SetData(oechem.OEGetTag('structure'), struct_out.encode(structure))
            oechem.OEWriteConstMolecule(ofs, mol)

        sys_in = OpenMMSystemInput('sys_in')
        struct_in = ParmEdStructureInput('struct_in')
        newmol = oechem.OEMol()
        with oechem.oemolistream('empty-mol.oeb.gz') as ifs:
            oechem.OEReadMolecule(ifs,newmol)
            newstruct = struct_in.decode(newmol.GetData(oechem.OEGetTag('structure')))
            newsys = sys_in.decode(newmol.GetData(oechem.OEGetTag('system')))
            self.assertEqual(newsys.getNumParticles(), 0)
            self.assertEqual(len(newstruct.atoms), 0)
        os.remove('empty-mol.oeb.gz')

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
        self.cube.args.steps = 1000
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def tearDown(self):
        self.runner.finalize()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        mol = oechem.OEMol()
        fname = get_data_filename('9PC1X-complex.oeb.gz')
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
        self.assertTrue(outmol.HasData(oechem.OEGetTag('structure')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('system')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('state')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('log')))

    def test_failure(self):
        pass

if __name__ == "__main__":
        unittest.main()
