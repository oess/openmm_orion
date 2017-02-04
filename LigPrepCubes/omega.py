from __future__ import print_function, division, unicode_literals
"""
Copyright (C) 2014, 2015 OpenEye Scientific Software
"""
import random

from openeye import oechem
from openeye import oeomega

from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube, MoleculeInputPort,
    StringParameter, MoleculeOutputPort, BinaryOutputPort
)
from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort


class OEOmegaConfGen(OEMolComputeCube):

    classification = [["OpenEye", "Conformer Generation"]]
    tags = "OMEGA,Conformer Generation".split(",")
    title = "OMEGA"

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    maxConfs = parameter.IntegerParameter('maxConfs', default=200,
                                          help_text="MaxConfs parameter for Omega")
    eWindow = parameter.IntegerParameter('eWindow', default=10,
                                         help_text="eWindow parameter for Omega")
    rms = parameter.DecimalParameter('rms', default=0.5,
                                     help_text="RMS Threshold parameter for Omega")
    maxSearchTime = parameter.IntegerParameter('maxSearchTime', default=120,
                                               help_text="Max Search Time parameter for Omega")
    use_ROC = parameter.BooleanParameter('use_ROC', default=True,
                                         help_text='Output molecules with Rotor-offset-compressed coords.')

    def begin(self):
        self.omega = oeomega.OEOmega()
        self.omega.SetMaxConfs(int(self.args.maxConfs))
        self.omega.SetEnergyWindow(float(self.args.eWindow))
        self.omega.SetRMSThreshold(float(self.args.rms))
        self.omega.SetMaxSearchTime(float(self.args.maxSearchTime))
        self.omega.SetRotorOffset(self.args.use_ROC)

    def process(self, mol, port):
        copy_mol = oechem.OEMol(mol)
        if self.omega(mol):
            if not mol.IsValid():
                self.log.error('omega returned empty molecule: {}'.format(copy_mol.GetTitle()))
                self.failure.emit(copy_mol)
            else:
                nconfs = mol.GetMaxConfIdx()
                self.log.info('OMEGA generated {} conformers for {}'.format(nconfs, mol.GetTitle()))
                self.success.emit(mol)
        else:
            self.failure.emit(copy_mol)


def _gen_flipper_wart(mol):
    res = ""
    for atom in mol.GetAtoms():
        if atom.HasData("__undef__"):
            cip = oechem.OEPerceiveCIPStereo(mol, atom)
            if cip == oechem.OECIPAtomStereo_R:
                res += "R"
            elif cip == oechem.OECIPAtomStereo_S:
                res += "S"
            atom.DeleteData("__undef__")
    for bond in mol.GetBonds():
        if bond.HasData("__undef__"):
            cip = oechem.OEPerceiveCIPStereo(mol, bond)
            if cip == oechem.OECIPBondStereo_E:
                res += "E"
            elif cip == oechem.OECIPBondStereo_Z:
                res += "Z"
            bond.DeleteData("__undef__")
    return res


def _label_unspecified(mol):
    for atom in mol.GetAtoms():
        if atom.IsChiral() and not atom.HasStereoSpecified():
            atom.SetBoolData("__undef__", True)
    for bond in mol.GetBonds():
        if bond.IsChiral() and not bond.HasStereoSpecified():
            bond.SetBoolData("__undef__", True)


class OEFlipStereo(OEMolComputeCube):

    title = 'Enumerate Stereo'

    classification = [["OpenEye", "Conformer Generation"]]

    add_label = parameter.BooleanParameter(
        'add_label',
        default=True,
        help_text='Add an annotation the title to show the enumerated stereo')
    max_centers = parameter.IntegerParameter(
        'max_centers',
        default=12,
        help_text='Max number of undefined stereo atoms/bonds to enumerate')

    def process(self, mol, port):
        if self.args.add_label:
            _label_unspecified(mol)
        for new_mol in oeomega.OEFlipper(mol):
            if self.args.add_label:
                wart = _gen_flipper_wart(new_mol)
                if len(wart) > 0:
                    new_mol.SetTitle(new_mol.GetTitle() + "_" + wart)
            res = oechem.OEMol(new_mol)
            # TODO: make sure this works on single conf molecules
            oechem.OECopySDData(res.GetActive(), mol.GetActive())
            self.success.emit(res)


class OEMakeMol3D(ParallelOEMolComputeCube):
    title = 'Generate Single Conformer'
    classification = [['Conformer Generation']]

    def begin(self):
        self.omega = oeomega.OEOmega()
        self.omega.SetMaxConfs(1)

    def process(self, mol, port):
        if self.omega(mol):
            self.success.emit(mol)
        else:
            self.failure.emit(mol)


class ConfSplitter(ParallelOEMolComputeCube):

    title = 'Multiplex Molecules with 10, 50 and 200 Conformers'
    classification = [["Conformer Generation"]]

    two_hundred_confs = MoleculeOutputPort('two_hundred_confs')
    fifty_confs = MoleculeOutputPort('fifty_confs')
    ten_confs = MoleculeOutputPort('ten_confs')

    def trim_confs(self, input_mol, num_confs):
        mol = oechem.OEMol(input_mol)
        for idx, conf in enumerate(mol.GetConfs()):
            if idx >= num_confs:
                mol.DeleteConf(conf)
        return mol

    def process(self, mol, port):
        mol200 = self.trim_confs(mol, 200)
        self.two_hundred_confs.emit(mol200)
        mol50 = self.trim_confs(mol200, 50)
        self.fifty_confs.emit(mol50)
        mol10 = self.trim_confs(mol50, 10)
        self.ten_confs.emit(mol10)


class OEOrionDBOStreamCube(SinkCube):
    """
    A sink cube that writes molecules to a file
    """
    title = 'Multiplex Molecules to files with 10, 50 and 200 Conformers'
    classification = [["Conformer Generation"], ["Output"]]
    intake = MoleculeInputPort("intake")
    prefix = StringParameter('prefix', required=True,
                             help_text='Prefix for the set out output files.')

    def begin(self):
        self.main_file = '{}_maxconfs200.oeb.gz'.format(self.args.prefix)
        self.ofs200 = oechem.oemolostream(self.main_file)
        self.ofs50 = oechem.oemolostream('{}_maxconfs50.oeb.gz'.format(self.args.prefix))
        self.ofs10 = oechem.oemolostream('{}_maxconfs10.oeb.gz'.format(self.args.prefix))

    def trim_confs(self, input_mol, num_confs):
        mol = oechem.OEMol(input_mol)
        for idx, conf in enumerate(mol.GetConfs()):
            if idx >= num_confs:
                mol.DeleteConf(conf)
        return mol

    def write(self, mol, port):
        mol200 = self.trim_confs(mol, 200)
        oechem.OEWriteMolecule(self.ofs200, mol200)
        mol50 = self.trim_confs(mol200, 50)
        oechem.OEWriteMolecule(self.ofs50, mol50)
        mol10 = self.trim_confs(mol50, 10)
        oechem.OEWriteMolecule(self.ofs10, mol10)

    def end(self):
        self.ofs200.close()
        self.ofs50.close()
        self.ofs10.close()

        # create 50K random set
        randomofs = oechem.oemolostream('{}_50Ksubset.oeb.gz'.format(self.args.prefix))
        moldb = oechem.OEMolDatabase(self.main_file)
        indices = range(moldb.NumMols())
        random.shuffle(indices)
        for index in indices[:50000]:
            moldb.GetMolecule(randomofs, index)
        randomofs.close()
