from __future__ import unicode_literals
from floe.api import OutputPort, InputPort
from floe.constants import BYTES
from simtk import openmm, unit
import io, base64, parmed
try:
    import cPickle as pickle
except ImportError:
    import pickle
