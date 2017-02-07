from __future__ import unicode_literals
from floe.api import OutputPort, InputPort
from floe.constants import BYTES
from simtk import openmm, unit
import io, base64, parmed
try:
    import cPickle as pickle
except ImportError:
    import pickle


class OpenMMSystemSerializationMixin(object):

    def encode(self, openmm_system):
        return openmm.XmlSerializer.serialize(openmm_system).encode()

    def decode(self, serialized_system):
        return openmm.XmlSerializer.deserialize(serialized_system)#.decode()

class OpenMMSystemInput(OpenMMSystemSerializationMixin, InputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")
class OpenMMSystemOutput(OpenMMSystemSerializationMixin, OutputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")

class ParmEdStructureSerializationMixin(object):

    def encode(self, structure):
        pkl_dict = pickle.dumps(structure.__getstate__())
        return base64.b64encode(pkl_dict)

    def decode(self, encoded_structure):
        decoded_structure = base64.b64decode(encoded_structure)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        return struct

class ParmEdStructureInput(ParmEdStructureSerializationMixin, InputPort):
    PORT_TYPES = (BYTES, "ParmEdStructure")
class ParmEdStructureOutput(ParmEdStructureSerializationMixin, OutputPort):
    PORT_TYPES = (BYTES, "ParmEdStructure")
