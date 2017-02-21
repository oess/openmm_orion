# -*- coding: utf-8 -*-
import io, os, parmed, mdtraj
from sys import stdout
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
from simtk.openmm import app
from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )
import io, base64, parmed
try:
    import cPickle as pickle
except ImportError:
    import pickle
# Prevents repeated downloads of the same Dataset
download_cache = {}


class OEPackMol(object):
    """A class designated to handle the packing/unpacking Python objects to the
    OEMol as generic data."""

    #def __init__(self, molecule):
    #    self.molecule = molecule
        #self.tag_data = dict()

    def getTags(molecule):
        return list(molecule.GetData().keys())

    def getData(molecule, tag):
        return molecule.GetData(oechem.OEGetTag(str(tag)))

    def encodeOpenMM(mm_obj):
        #Supports State, System, Integrator, Forcefield
        return openmm.XmlSerializer.serialize(mm_obj).encode()

    def decodeOpenMM(data):
        return openmm.XmlSerializer.deserialize(data)

    def encodeStruct(structure):
        pkl_dict = pickle.dumps(structure.__getstate__())
        return base64.b64encode(pkl_dict)

    def decodeStruct(data):
        decoded_structure = base64.b64decode(data)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        return struct

    def getState(simulation):


        return state

    def encodeSimData(idtag, simulation):
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

        trajfname = idtag+'-simulation.nc'
        traj = mdtraj.load(trajfname, top=mdtraj.Topology.from_openmm(topology))

        structure = parmed.openmm.load_topology(topology, system,
                                   xyz=state.getPositions())

        tag_data['state'] = OEPackMol.encodeOpenMM(state)
        tag_data['pdb'] = OEPackMol.encodePDBFile(state, topology)
        tag_data['traj'] = OEPackMol.encodeTrajectory(traj)
        tag_data['log'] = open(idtag+'-simulation.log').read()
        tag_data['structure'] = OEPackMol.encodeStruct(structure)
        return tag_data

    def encodeTrajectory(trajectory):
        pkl_traj = pickle.dumps(trajectory)
        return base64.b64encode(pkl_traj)

    def decodeTrajectory(data):
        decoded_traj = base64.b64decode(data)
        return pickle.loads(decoded_traj)

    def encodePDBFile(state, topology):
        f = io.StringIO()
        app.PDBFile.writeFile(topology, state.getPositions(), f)
        return f.getvalue()

    def decodePDBFile(data):
        return app.PDBFile(io.StringIO(data))

    @staticmethod
    def checkTags(molecule, tags):
        oetags = OEPackMol.getTags(molecule)
        intersect = list( set(tags) & set(oetags) )
        diff = list( set(tags) - set(oetags) )
        if diff:
            raise RuntimeError('Missing {} in tagged data'.format(diff))
        else:
            #print('Found tags: {}'.format(intersect))
            return True

    @classmethod
    def unpack(cls, molecule):
        tag_data = {}
        tags = cls.getTags(molecule)
        for tag in tags:
            data = cls.getData(molecule, tag)
            if 'state' == tag:
                data = cls.decodeOpenMM(data)
            if 'structure' == tag:
                data = cls.decodeStruct(data)
            if 'pdb' == tag:
                data = cls.decodePDBFile(data)
            if 'traj' == tag:
                data = cls.decodeTrajectory(data)
            tag_data[tag] = data
        return tag_data

    @classmethod
    def pack(cls, molecule, data):
        try:
            idtag = cls.getData(molecule, 'idtag')
        except Exception as e:
            print(e)
            pass

        if isinstance(data, parmed.structure.Structure):
            molecule.SetData(oechem.OEGetTag('structure'), cls.encodeStruct(data))

        if isinstance(data, openmm.app.simulation.Simulation):
            tag_data = cls.encodeSimData(idtag, data)
            for k,v in tag_data.items():
                molecule.SetData(oechem.OEGetTag(k), v)

        return molecule

def cleanup(tmpfiles):
    for tmp in tmpfiles:
        try:
            os.remove(tmp)
        except Exception as e:
            pass

def setReporters(reportInterval, totalSteps, outfname):
    progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                        reportInterval=reportInterval, totalSteps=totalSteps,
                                        time=True, speed=True, progress=True,
                                        elapsedTime=True, remainingTime=True)

    state_reporter = app.StateDataReporter(outfname+'.log', separator="\t",
                                        reportInterval=reportInterval,
                                        step=True,
                                        potentialEnergy=True, totalEnergy=True,
                                        volume=True, temperature=True)

    #traj_reporter = mdtraj.reporters.HDF5Reporter(self.outfname+'.h5', self.args.reporter_interval)
    traj_reporter = mdtraj.reporters.NetCDFReporter(outfname+'.nc', reportInterval)
    ##dcd_reporter = app.dcdreporter.DCDReporter(self.outfname+'.dcd', self.args.reporter_interval)
    reporters = [state_reporter, progress_reporter, traj_reporter]
    return reporters

def genSimFromStruct(structure, temperature):
    system = structure.createSystem(nonbondedMethod=app.PME,
                                    nonbondedCutoff=10.0*unit.angstroms,
                                    constraints=app.HBonds)

    integrator = openmm.LangevinIntegrator(temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
    simulation = app.Simulation(structure.topology, system, integrator)
    simulation.context.setPositions(structure.positions)
    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    platform = simulation.context.getPlatform().getName()
    print('Running OpenMMSimulation on Platform {}'.format(platform))
    return simulation

def minimizeSimulation(simulation):
        init = simulation.context.getState(getEnergy=True)
        print('Initial energy is {}'.format(init.getPotentialEnergy()))
        print('Minimizing system...')
        #f = io.StringIO()
        #rep = parmed.openmm.EnergyMinimizerReporter(f, volume=True)
        simulation.minimizeEnergy()
        #rep.report(simulation)
        #print(f.getvalue())
        st = simulation.context.getState(getPositions=True,getEnergy=True)
        print('Minimized energy is {}'.format(st.getPotentialEnergy()))
        return simulation

        #minfname = 'output/{}-minimized.pdb'.format(idtag)
        #with open(minfname, 'w') as minout:
        #    app.PDBFile.writeFile(simulation.topology, min_state.getPositions(), minout)

        #with open(minfname, 'rb') as minin:
        #    outmol.SetData(oechem.OEGetTag('pdb') , minin.read())
        #os.remove(minfname)


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
