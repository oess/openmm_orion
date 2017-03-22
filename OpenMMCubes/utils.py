import io, os, base64, parmed, mdtraj, pdbfixer, glob, tarfile
import numpy as np
from sys import stdout
from tempfile import NamedTemporaryFile
from openeye import oechem
from floe.api.orion import in_orion, StreamingDataset, upload_file
from simtk import unit, openmm
from simtk.openmm import app
try:
    import cPickle as pickle
except ImportError:
    import pickle
# Prevents repeated downloads of the same Dataset
download_cache = {}

class PackageOEMol(object):
    """
    Class designated to handle the packing/unpacking of the ParmEd Structure,
    OpenMM State, and the log file for the state reporter. Objects are attached
    to the OEMol as generic data. Attachment of the ParmEd Structure is required
    in order to preserve the SMIRFF parameters during complex generation.
    """

    def getTags(molecule):
        return list(molecule.GetData().keys())

    def getData(molecule, tag):
        return molecule.GetData(oechem.OEGetTag(str(tag)))

    def updateSimIdx(molecule):
        """Place a counter on the OEMol to track the number of restarts"""
        if molecule.HasData(oechem.OEGetTag('SimIdx')):
            idx = molecule.GetData(oechem.OEGetTag('SimIdx'))
            ndx = int(idx)+1
            if  ndx <= 9:
                simidx = '0%s' % ndx
            else:
                simidx = str(ndx)
        else:
            simidx = '01'
        molecule.SetData(oechem.OEGetTag('SimIdx'), simidx)
        return molecule

    def encodeOpenMM(mm_obj):
        """Serialize the OpenMM Objects to BYTES.
        Supports State, System, Integrator, Forcefield objects.
        """
        return openmm.XmlSerializer.serialize(mm_obj).encode()

    def decodeOpenMM(data):
        """Deserialize the bytes reperesentation of the OpenMM Objects to STR"""
        return openmm.XmlSerializer.deserialize(data)

    def encodePyObj(py_obj):
        """Encode Python Objects/ParmEd Structure to BYTES (base64)"""
        pkl_obj = pickle.dumps(py_obj)
        return base64.b64encode(pkl_obj)

    def decodePyObj(data):
        """Decode the Base64 encoded Python Object/ParmEd Structure"""
        decoded_obj = base64.b64decode(data)
        return pickle.loads(decoded_obj)

    def encodeStruct(structure):
        """Encode the ParmEd Structure by retrieving the current state as a dict.
        Dict contains objects and cannot be encoded directly by pickle."""
        struct_dict = structure.__getstate__()
        return PackageOEMol.encodePyObj(struct_dict)

    def decodeStruct(data):
        """Decode the Base64 encoded Structure dict, then
        load into an empty Structure object."""
        struct = parmed.structure.Structure()
        struct.__setstate__(PackageOEMol.decodePyObj(data))
        return struct

    def encodeSimData(simulation):
        """Pulls the OpenMM State object and log file reproting energies from the
        OpenMM Simulation object. Generates a new ParmEd Structure from the
        State of the simulation."""
        tag_data = {}
        topology = simulation.topology
        system = simulation.context.getSystem()
        state = simulation.context.getState(getPositions=True,
                                            getVelocities=True,
                                            getParameters=True,
                                            getForces=True,
                                            getParameterDerivatives=True,
                                            getEnergy=True,
                                            enforcePeriodicBox=True)

        # Get the proper log file name from the reporters, not stdout
        for rep in simulation.reporters:
            if isinstance(rep, app.statedatareporter.StateDataReporter):
                if rep._out is stdout:
                    pass
                else:
                    logfname = rep._out.name

        # Generate the ParmEd Structure
        structure = parmed.openmm.load_topology(topology, system,
                                   xyz=state.getPositions())

        # Return dictionary with encoded data
        tag_data['State'] = PackageOEMol.encodeOpenMM(state)
        tag_data['Structure'] = PackageOEMol.encodeStruct(structure)
        with open(logfname) as log:
            tag_data['Log'] = log.read()
        return tag_data

    def checkSDData(molecule):
        """ Returns a dictionary of the SD Data from the OEMol """
        sd_data = {}
        for dp in oechem.OEGetSDDataPairs(molecule):
            sd_data[dp.GetTag()] = dp.GetValue()
        return sd_data

    @staticmethod
    def checkTags(molecule, req_tags=None):
        """ Checks if OEMol has required data attached """
        oetags = PackageOEMol.getTags(molecule)
        intersect = list( set(req_tags) & set(oetags) )
        diff = list( set(req_tags) - set(oetags) )
        if diff:
            raise RuntimeError('Missing {} in tagged data'.format(diff))
        else:
            #print('Found tags: {}'.format(intersect))
            return True

    @classmethod
    def unpack(cls, molecule, tags=None):
        """ Decodes the data attached to OEMol. Return as dictionary """
        tag_data = {}
        # Default to decode all
        if tags is None:
            tags = cls.getTags(molecule)
        for tag in tags:
            data = cls.getData(molecule, tag)
            if 'State' == tag:
                data = cls.decodeOpenMM(data)
            if 'Structure' == tag:
                data = cls.decodeStruct(data)
            if 'Log' == tag:
                data = io.StringIO(data)
            tag_data[tag] = data
        return tag_data

    @classmethod
    def dump(cls, molecule, tags=None, outfname=None, tarxz=True):
        """ Writes the data attached to OEMol to files on disc.
        Create a tar archive (xz-compressed) of files."""

        tag_data = {}
        totar = []
        if tags is None:
            tag_data = cls.unpack(molecule)
        else:
            tag_data = cls.unpack(molecule, tags=tags)
        if outfname is None:
            raise Exception('Require an output file name.')

        # Dump the SD Data
        sd_data = cls.checkSDData(molecule)
        sdtxt = outfname+'-sd.txt'
        with open(sdtxt, 'w') as f:
            for k,v in sd_data.items():
                f.write('{} : {}\n'.format(k,v))
        totar.append(sdtxt)

        print('Dumping data from: %s' % outfname)
        for tag, data in tag_data.items():
            if isinstance(data, parmed.structure.Structure):
                pdbfname = outfname+'.pdb'
                print("\tStructure to %s" % pdbfname)
                data.save(pdbfname, overwrite=True)
                totar.append(pdbfname)
            if isinstance(data, openmm.openmm.State):
                statefname = outfname+'-state.xml'
                print('\tState to %s' % statefname)
                with open(statefname, 'w') as f:
                    f.write(openmm.XmlSerializer.serialize(data))
                if tarxz:
                    totar.append(statefname)
            if isinstance(data, io.StringIO):
                enelog = outfname+'.log'
                print('\tLog to %s' % enelog)
                with open(enelog, 'w') as f:
                    f.write(data.getvalue())
                if tarxz:
                    totar.append(enelog)
        if tarxz:
            tarname = outfname+'.tar.xz'
            print('Creating tarxz file: {}'.format(tarname))

            trajfname = outfname+'.nc'
            if os.path.isfile(trajfname):
               totar.append(trajfname)
               print('Adding {} to {}'.format(trajfname, tarname))
            else:
               print('Could not find {}'.format(trajfname))

            tar = tarfile.open(tarname, "w:xz")
            for name in totar:
                tar.add(name)
            tar.close()

            #### MUST upload tar file directly back to Orion or they disappear.
            upload_file(tarname, tarname, tags=['TAR'])
            # Clean up files that have been added to tar.
            cleanup(totar)

    @classmethod
    def pack(cls, molecule, data):
        """ Encodes the ParmEd Structure or if provided the OpenMM Simulation object,
        this will extract the State and the log file from the state reporter and
        attach them to the OEMol as generic data. Returns the OEMol with attached data."""

        tag_data = {}
        # Attach (base64) encoded ParmEd Structure.
        if isinstance(data, parmed.structure.Structure):
            molecule.SetData(oechem.OEGetTag('Structure'), cls.encodeStruct(data))

        # Attach the encoded OpenMM State and log file from the Simulation.
        if isinstance(data, openmm.app.simulation.Simulation):
            tag_data = cls.encodeSimData(data)
            for k,v in tag_data.items():
                molecule.SetData(oechem.OEGetTag(k), v)
        return molecule

def cleanup(tmpfiles):
    for tmp in tmpfiles:
        try:
            os.remove(tmp)
        except Exception as e:
            pass

def get_data_filename(package_root, relative_path):
    """Get the full path of the files installed in python packages or included
    in this package.

    Parameters
    ----------
    package_root : str (i.e examples or smirff99Frosst)
        Name of the included/installed python package
    relative_path: str (i.e. toluene.pdb or smirff99Frosst.ffxml)
        Path to the file within the python package

    Returns
    --------
    fn : str (i.e examples/data/toluene.pdb or smirff99Frosst/smirff99Frosst.ffxml)
        Full path to the file within the python package
    """

    from pkg_resources import resource_filename
    fn = resource_filename(package_root, os.path.join(relative_path))
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
