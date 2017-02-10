# -*- coding: utf-8 -*-
import os

from tempfile import NamedTemporaryFile

from openeye.oechem import(
    oemolostream, OEWriteConstMolecule
)
from openeye import oechem
from openeye.oedocking import OEWriteReceptorFile
import numpy as np
from floe.api.orion import in_orion, StreamingDataset
from simtk.openmm.app import Topology
from simtk.openmm.app.element import Element
from simtk import unit, openmm
# Prevents repeated downloads of the same Dataset
download_cache = {}

def get_data_filename(relative_path):
    """Get the full path to one of the reference files in testsystems.
    In the source distribution, these files are in ``examples/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.
    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).
    """

    from pkg_resources import resource_filename
    fn = resource_filename('examples', os.path.join('data', relative_path))
    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)
    return fn

def getPositionsFromOEMol(molecule):
    positions = unit.Quantity(
        np.zeros([molecule.NumAtoms(), 3], np.float32), unit.angstroms)
    coords = molecule.GetCoords()
    for index in range(molecule.NumAtoms()):
        positions[index, :] = unit.Quantity(coords[index], unit.angstroms)
    return positions

def combinePostions(proteinPositions, molPositions):
    # Concatenate positions arrays (ensures same units)
    positions_unit = unit.angstroms
    positions0_dimensionless = np.array(proteinPositions / positions_unit)
    positions1_dimensionless = np.array(molPositions / positions_unit)
    coordinates = np.vstack(
        (positions0_dimensionless, positions1_dimensionless))
    natoms = len(coordinates)
    positions = np.zeros([natoms, 3], np.float32)
    for index in range(natoms):
            (x, y, z) = coordinates[index]
            positions[index, 0] = x
            positions[index, 1] = y
            positions[index, 2] = z
    positions = unit.Quantity(positions, positions_unit)
    return positions

def download_dataset_to_file(dataset_id):
    """
    Used to retrieve a dataset either from Orion or from the local machine
    """
    if in_orion():
        if dataset_id in download_cache:
            return download_cache[dataset_id]
        if os.path.isfile(dataset_id):
            download_cache[dataset_id] = dataset_id
            return dataset_id
        tmp = NamedTemporaryFile(suffix=".oeb.gz", delete=False)
        stream = StreamingDataset(dataset_id, input_format=".oeb.gz")
        stream.download_to_file(tmp.name)
        download_cache[dataset_id] = tmp.name
        return tmp.name
    else:
        return dataset_id


def dump_query(prefix, name, qmol, receptor):
    """
    Writes the Molecule or receptor out to file on the machine
    """
    tag = "{0}_{1}.query".format(prefix, name)
    query_file = "{0}.oeb.gz".format(tag)
    with oemolostream(query_file) as ofs:
        OEWriteConstMolecule(ofs, qmol)
    if receptor.IsValid():
        receptor_file = "{0}.receptor.oeb.gz".format(tag)
        OEWriteReceptorFile(receptor, receptor_file)
    return tag, query_file
