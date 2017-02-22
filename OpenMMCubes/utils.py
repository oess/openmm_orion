import io, os, base64, parmed, mdtraj, pdbfixer, glob
import numpy as np
from sys import stdout
from tempfile import NamedTemporaryFile
from openeye import oechem
from floe.api.orion import in_orion, StreamingDataset
from simtk import unit, openmm
from simtk.openmm import app
from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
try:
    import cPickle as pickle
except ImportError:
    import pickle
# Prevents repeated downloads of the same Dataset
download_cache = {}


class PackageOEMol(object):
    """A class designated to handle the packing/unpacking Python objects to the
    OEMol as generic data."""

    def getTags(molecule):
        return list(molecule.GetData().keys())

    def getData(molecule, tag):
        return molecule.GetData(oechem.OEGetTag(str(tag)))

    def encodeOpenMM(mm_obj):
        #Supports State, System, Integrator, Forcefield
        return openmm.XmlSerializer.serialize(mm_obj).encode()

    def decodeOpenMM(data):
        return openmm.XmlSerializer.deserialize(data)

    def encodePyObj(py_obj):
        pkl_obj = pickle.dumps(py_obj)
        return base64.b64encode(pkl_obj)

    def decodePyObj(data):
        decoded_obj = base64.b64decode(data)
        return pickle.loads(decoded_obj)

    def encodeStruct(structure):
        struct_dict = structure.__getstate__()
        return PackageOEMol.encodePyObj(struct_dict)

    def decodeStruct(data):
        struct = parmed.structure.Structure()
        struct.__setstate__(PackageOEMol.decodePyObj(data))
        return struct

    def encodeSimData(simulation, **opt):
        tag_data = {}
        #Close open files for reporters.
        for rep in simulation.reporters:
            try:
                rep.close()
            except:
                pass

        topology = simulation.topology
        system = simulation.context.getSystem()
        state = simulation.context.getState(getPositions=True,
                                            getVelocities=True,
                                            getParameters=True,
                                            getForces=True,
                                            getParameterDerivatives=True,
                                            getEnergy=True,
                                            enforcePeriodicBox=True)

        traj = mdtraj.load(opt['outfname']+'.h5')
        structure = parmed.openmm.load_topology(topology, system,
                                   xyz=state.getPositions())

        tag_data['state'] = PackageOEMol.encodeOpenMM(state)
        tag_data['traj'] = PackageOEMol.encodePyObj(traj)
        tag_data['log'] = open(opt['outfname']+'.log').read()
        tag_data['structure'] = PackageOEMol.encodeStruct(structure)

        tmpfiles = [ opt['outfname']+'.log', opt['outfname']+'.h5' ]
        cleanup(tmpfiles)
        return tag_data

    @staticmethod
    def checkTags(molecule, tags):
        oetags = PackageOEMol.getTags(molecule)
        intersect = list( set(tags) & set(oetags) )
        diff = list( set(tags) - set(oetags) )
        if diff:
            raise RuntimeError('Missing {} in tagged data'.format(diff))
        else:
            #print('Found tags: {}'.format(intersect))
            return True

    @classmethod
    def unpack(cls, molecule, tags=None):
        tag_data = {}
        if tags is None:
            tags = cls.getTags(molecule)
        for tag in tags:
            data = cls.getData(molecule, tag)
            if 'state' == tag:
                data = cls.decodeOpenMM(data)
            if 'structure' == tag:
                data = cls.decodeStruct(data)
            if 'traj' == tag:
                data = cls.decodePyObj(data)
            if 'log' == tag:
                data = io.StringIO(data)
            tag_data[tag] = data
        return tag_data

    @classmethod
    def dump(cls, molecule, tags=None, **opt):
        tag_data = {}
        log = opt['logger']
        traj = { 'HDF5' : opt['outfname']+'.h5',
                 'DCD' : opt['outfname']+'.dcd',
                 'NetCDF' : opt['outfname']+'.nc' }

        if tags is None:
            tag_data = cls.unpack(molecule)
        else:
            tag_data = cls.unpack(molecule, tags=tags)

        log.info('Dumping %s data...' % opt['outfname'])
        for tag, data in tag_data.items():
            if isinstance(data, parmed.structure.Structure):
                pdbfname = opt['outfname']+'.pdb'
                log.info("\tStructure to %s" % pdbfname)
                data.save(pdbfname, overwrite=True)

            if isinstance(data, mdtraj.core.trajectory.Trajectory):
                trajfname = traj.get(opt['trajectory_filetype'])
                log.info("\tTrajectory to %s" % trajfname)
                data.save(trajfname)

            if isinstance(data, openmm.openmm.State):
                statefname = opt['outfname']+'-state.xml'
                log.info('\tState to %s' % statefname)
                with open(statefname, 'w') as f:
                    f.write(openmm.XmlSerializer.serialize(data))

            if isinstance(data, io.StringIO):
                enelog = opt['outfname']+'.log'
                log.info('\tLog to %s' % enelog)
                with open(enelog, 'w') as f:
                    f.write(data.getvalue())

    @classmethod
    def pack(cls, molecule, data, tags=None, dump=False, **opt):
        if isinstance(data, parmed.structure.Structure):
            molecule.SetData(oechem.OEGetTag('structure'), cls.encodeStruct(data))

        if isinstance(data, openmm.app.simulation.Simulation):
            tag_data = cls.encodeSimData(data,**opt)
            for k,v in tag_data.items():
                molecule.SetData(oechem.OEGetTag(k), v)

        if dump:
            if tags:
                cls.dump(molecule, tags, **opt)
            else:
                cls.dump(molecule, **opt)

        return molecule

def cleanup(tmpfiles):
    for tmp in tmpfiles:
        try:
            os.remove(tmp)
        except Exception as e:
            pass

def getPositionsFromOEMol(molecule):
    positions = unit.Quantity(
        np.zeros([molecule.NumAtoms(), 3], np.float32), unit.angstroms)
    coords = molecule.GetCoords()
    for index in range(molecule.NumAtoms()):
        positions[index, :] = unit.Quantity(coords[index], unit.angstroms)
    return positions

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
