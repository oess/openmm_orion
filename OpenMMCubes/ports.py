from __future__ import unicode_literals
from floe.api import OutputPort, InputPort
from floe.constants import BYTES
from simtk import openmm, unit
from lzma import compress, decompress



class OpenMMSystemSerializationMixin(object):

    def encode(self, openmm_system):
        #return compress( openmm.XmlSerializer.serialize(openmm_system).encode() )
        return openmm.XmlSerializer.serialize(openmm_system).encode()

    def decode(self, serialized_system):
        #return openmm.XmlSerializer.deserialize( decompress(serialized_system).decode() )
        return openmm.XmlSerializer.deserialize(serialized_system)#.decode()


class OpenMMSystemInput(OpenMMSystemSerializationMixin, InputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")

class OpenMMSystemOutput(OpenMMSystemSerializationMixin, OutputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")
