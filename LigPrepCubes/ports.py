from __future__ import unicode_literals
from floe.api import OutputPort, InputPort, Port
from floe.constants import BYTES
from openeye.oechem import (
    oemolistream, OEMol, oemolostream, OEWriteMolecule, oeosstream,
    OEWriteMolReturnCode_Success, OEFormat_OEB,
    OEReadMolecule, OEGetTag)
#from floe.constants import BYTES
import io
try:
    import cPickle as pickle
except ImportError:
    import pickle

class MoleculeSerializerMixin(object):

    def __init__(self, *args, **kwargs):
        super(MoleculeSerializerMixin, self).__init__(*args, **kwargs)
        self._ifs = oemolistream()
        self._ifs.SetFormat(OEFormat_OEB)
        self._ifs.Setgz(True)
        errs = oeosstream()
        self._ofs = oemolostream(errs, False)
        self._ofs.openstring()
        self._ofs.Setgz(True)
        self._ofs.SetFormat(OEFormat_OEB)

    def encode(self, *mols):
        """
        By default, serializes molecules as gzipped oeb files in a string
        """
        for mol in mols:
            code = OEWriteMolecule(self._ofs, mol)
            if code != OEWriteMolReturnCode_Success:
                raise RuntimeError("Unable to encode mol: {}".format(code))
        res = self._ofs.GetString()
        self._ofs.close()
        self._ofs.openstring()
        return res

    def decode(self, mol_data):
        """
        By default, deserializes data into molecules for use in the cube
        """
        mol = OEMol()
        if type(mol_data) == OEMol:
            return mol_data
        if not self._ifs.openstring(mol_data):
            raise RuntimeError("Failed to open string")
        if not OEReadMolecule(self._ifs, mol):
            print("Unable to decode molecule")
        self._ifs.close()
        return mol


class MoleculePortSerializer(MoleculeSerializerMixin, Port):

    OEB_GZ = '.oeb.gz'
    PORT_TYPES = (BYTES, OEB_GZ)
    FORMAT = OEB_GZ


class CustomMoleculeInputPort(InputPort, MoleculePortSerializer):
    pass


class CustomMoleculeOutputPort(OutputPort, MoleculePortSerializer):
    pass
