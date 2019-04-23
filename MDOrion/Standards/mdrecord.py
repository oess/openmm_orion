from MDOrion.Standards import (Fields,
                               MDFileNames,
                               MDEngines,
                               MDStageTypes)

from MDOrion.Standards import utils

import parmed

import copy

import os

import tarfile

import tempfile

from tempfile import TemporaryDirectory

from datarecord import (Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from openeye import oechem

import pickle

import shutil

from orionclient.types import (ShardCollection,
                               Shard)

from orionclient.session import (in_orion,
                                 APISession)

import glob


def mdstages(f):

    def wrapper(*pos, **named):
        mdrec = pos[0]

        if not mdrec.rec.has_field(Fields.md_stages):
            raise ValueError("The MD record does not have MD stages")

        return f(*pos, **named)

    return wrapper


def stage_system(f):
    def wrapper(*pos, **named):

        mdrec = pos[0]

        if 'stg_name' not in named.keys():
            named['stg_name'] = 'last'

        stg_name = named['stg_name']

        stage = mdrec.get_stage_by_name(stg_name)

        stage_name = stage.get_value(Fields.stage_name)

        dir_stage = mdrec.processed[stage_name]

        if not dir_stage:

            dir_stage = tempfile.mkdtemp(prefix=stage_name + '_', dir=mdrec.cwd)

            file_id = stage.get_value(Fields.mddata)

            fn = utils.download_data(file_id, dir_stage, collection_id=mdrec.collection_id)

            with tarfile.open(fn) as tar:
                tar.extractall(path=dir_stage)

            mdrec.processed[stage_name] = dir_stage

        return f(*pos, **named)

    return wrapper


class MDDataRecord(object):

    def __init__(self, record, inplace=True):
        if inplace:
            self.rec = record
        else:
            self.rec = copy.deepcopy(record)

        if not self.rec.has_field(Fields.md_stages):
            self.processed = {}
        else:
            stages = self.rec.get_value(Fields.md_stages)
            self.processed = {stg.get_value(Fields.stage_name): False for stg in stages}

        if in_orion():
            if self.rec.has_field(Fields.collection):
                self.collection_id = self.rec.get_value(Fields.collection)

        else:
            self.collection_id = None

        self.cwd = tempfile.mkdtemp()

    def __del__(self):
        try:
            shutil.rmtree(self.cwd, ignore_errors=True)
        except OSError as e:
            print("Error: {} - {}".format(e.filename, e.strerror))

    @property
    def get_record(self):
        """
        This method returns the record

        Return:
        -------
        record: OERecord
            The record
        """
        return self.rec

    @property
    def get_primary(self):
        """
        This method returns the primary molecule present on the record

        Return:
        -------
        record: OEMol
            The Primary Molecule
        """

        if not self.rec.has_field(Fields.primary_molecule):
            raise ValueError("The Primary Molecule has not been found on the record")

        return self.rec.get_value(Fields.primary_molecule)

    def set_primary(self, primary_mol):
        """
        This method sets the primary molecule on the record

        Parameters:
        -----------
        primary_mol: OEMol
            The primary molecule to set on the record

        Return:
        -------
        record: Bool
            True if the primary molecule has been set on the record
        """

        if not isinstance(primary_mol, oechem.OEMol):
            raise ValueError("The Primary Molecule is not a valid OEMol: {}".format(primary_mol))

        self.rec.set_value(Fields.primary_molecule, primary_mol)

        return True

    @property
    def get_id(self):
        """
        This method returns the identification field ID present on the record

        Return:
        -------
            : Int
            The record ID
        """

        if not self.rec.has_field(Fields.id):
            raise ValueError("The ID Field has not been found on the record")

        return self.rec.get_value(Fields.id)

    @property
    def has_id(self):
        """
        This method checks if the identification field ID is present on the record

        Return:
        -------
            : Bool
            True if the ID field is resent on the record otherwise False
        """

        if not self.rec.has_field(Fields.id):
            return False
        else:
            return True

    def set_id(self, id):
        """
        This method sets the identification field ID on the record

        Parameters:
        -----------
        id: Int
            An identification integer for the record

        Return:
        -------
            : Bool
            True if the ID has been set on the record
        """

        if not isinstance(id, int):
            raise ValueError(" The id must be an integer: {}".format(id))

        self.rec.set_value(Fields.id, id)

        return True

    @property
    def get_title(self):
        """
        This method returns the system title present on the record

        Return:
        -------
            : String
            The system title tring
        """

        if not self.rec.has_field(Fields.title):
            raise ValueError("The Title Field has not been found on the record")

        return self.rec.get_value(Fields.title)

    @property
    def has_title(self):
        """
        This method checks if the Title field is present on the record

        Return:
        -------
            : Bool
            True if the Title field is resent on the record otherwise False
        """

        if not self.rec.has_field(Fields.title):
            return False
        else:
            return True

    def set_title(self, title):
        """
        This method sets the system Title field on the record

        Parameters:
        -----------
        title: String
            A string used to identify the  molecular system

        Return:
        -------
            : Bool
            True if the system Title has been set on the record
        """

        if not isinstance(title, str):
            raise ValueError(" The title must be a sting: {}".format(title))

        self.rec.set_value(Fields.title, title)

        return True

    def create_collection(self, name):
        """
        This method sets a collection field on the record to be used in Orion

        Parameters:
        -----------
        name: String
            A string used to identify in the Orion UI the record collection

        Return:
        -------
            : Bool
                True if the collection creation was successful
        """

        if in_orion():

            if self.rec.has_field(Fields.collection):
                raise ValueError("Collection field already present on the record")

            session = APISession

            collection = ShardCollection.create(session, name)

            self.rec.set_value(Fields.collection, collection.id)

            self.collection_id = collection.id

        else:
            return False

        return True

    @property
    @mdstages
    def get_last_stage(self):
        """
        This method returns the last MD stage of the MD record stages

        Return:
        -------
        record: OERecord
            The last stage of the MD record stages
        """

        return self.rec.get_value(Fields.md_stages)[-1]

    @mdstages
    def get_stage_by_name(self, stg_name='last'):
        """
        This method returns a MD stage selected by passing the string stage name. If the
        string "last" is passed (default) the last MD stage is returned. If multiple stages
        have the same name the first occurrence is returned. If no stage name has been found
        an exception is raised.

        Parameters:
        -----------
        stg_name: String
            The MD stage name

        Return:
        -------
        record: OERecord
            The MD stage selected by its name
        """
        stages = self.rec.get_value(Fields.md_stages)

        stg_names = []

        if stg_name == 'last':
            return stages[-1]

        for stage in stages:
            name = stage.get_value(Fields.stage_name)
            if name == stg_name:
                return stage
            else:
                stg_names.append(name)

        raise ValueError("The Stage name has not been found: {} available names: {}".format(stg_name, stg_names))

    @mdstages
    def get_stage_by_idx(self, idx):
        """
        This method returns a MD stage selected by passing the an index. If the index is not present
        an exception is raised.

        Parameters:
        -----------
        idx: Int
            The stage index to retrieve

        Return:
        -------
        record: OERecord
            The MD stage selected by its index
        """

        if idx > len(self.processed):
            raise ValueError("The selected stage index is greater than the md stages size {} > {}".
                             format(idx, len(self.processed)))

        return self.rec.get_value(Fields.md_stages)[idx]

    @mdstages
    def delete_stage_by_name(self, stg_name='last'):
        """
        This method deletes an MD stage selected by passing the string name. If the
        string "last" is passed (default) the last MD stage is deleted. If no stage name
        has been found an exception is raised.

        Parameters:
        -----------
        stg_name: String
            The MD stage name

        Return:
        -------
        record: Bool
            True if the deletion was successful
        """
        stages = self.rec.get_value(Fields.md_stages)

        stg_names = []

        if len(self.get_stages) == 1:

            stage = self.get_stage_by_idx(0)

            fid = stage.get_value(Fields.mddata)
            utils.delete_data(fid, collection_id=self.collection_id)

            if stage.get_value(Fields.trajectory) is not None:
                tid = stage.get_value(Fields.trajectory)
                utils.delete_file(tid)

            self.rec.delete_field(Fields.md_stages)
            self.processed = {}

            return True

        if stg_name == 'last':
            last_stage = stages[-1]
            name = last_stage.get_value(Fields.stage_name)
            fid = last_stage.get_value(Fields.mddata)
            utils.delete_data(fid, collection_id=self.collection_id)

            if last_stage.get_value(Fields.trajectory) is not None:
                tid = last_stage.get_value(Fields.trajectory)
                utils.delete_file(tid)

            del self.processed[name]
            del stages[-1]
            self.rec.set_value(Fields.md_stages, stages)

            return True
        else:
            for stage in stages:

                name = stage.get_value(Fields.stage_name)

                if name == stg_name:
                    fid = stage.get_value(Fields.mddata)
                    utils.delete_data(fid, collection_id=self.collection_id)

                    if stage.get_value(Fields.trajectory) is not None:
                        tid = stage.get_value(Fields.trajectory)
                        utils.delete_file(tid)

                    del self.processed[name]
                    stages.remove(stage)
                    self.rec.set_value(Fields.md_stages, stages)
                    return True
                else:
                    stg_names.append(name)

        raise ValueError("The Stage name has not been found: {} available names: {}".format(stg_name, stg_names))

    @mdstages
    def delete_stage_by_idx(self, idx):
        """
        This method deletes an MD stage selected by passing its index. If the stage index
        cannot be found an exception is raised.

        Parameters:
        -----------
        idx: Int
            The MD stage index

        Return:
        -------
        record: Bool
            True if the deletion was successful
        """

        stage = self.get_stage_by_idx(idx)
        name = stage.get_value(Fields.stage_name)

        return self.delete_stage_by_name(stg_name=name)

    @mdstages
    def has_stage_name(self, stg_name):
        """
        This method returns True if MD stage selected by passing the string name is present
        on the MD stage record otherwise False.

        Parameters:
        -----------
        name: String
            The MD stage name

        Return:
        -------
        record: Bool
            True if the MD stage name is present on the MD stages record otherwise False
        """

        stages = self.rec.get_value(Fields.md_stages)

        for stage in stages:
            name = stage.get_value(Fields.stage_name)
            if name == stg_name:
                return True

        return False

    @mdstages
    def get_stage_info(self, stg_name='last'):
        """
        This method returns the info related to the selected stage name. If no stage name is passed
        the last stage is selected.

        Parameters:
        -----------
        stg_name: String
            The MD stage name

        Return:
        -------
            : String
            The info associated with the selected MD stage
        """

        stage = self.get_stage_by_name(stg_name)

        return stage.get_value(Fields.log_data)

    @stage_system
    @mdstages
    def unpack_stage_system(self, stg_name='last'):
        pass

    @stage_system
    @mdstages
    def get_stage_state(self, stg_name='last'):
        """
        This method returns the MD State of the selected stage name. If no stage name is passed
        the last stage is selected

        Parameters:
        -----------
        stg_name: String
            The MD stage name
        Return:
        -------
            : OERecord
            The MD state of the selected MD stage
        """
        stage = self.get_stage_by_name(stg_name)

        stage_name = stage.get_value(Fields.stage_name)

        dir_stage = self.processed[stage_name]

        state_fn = os.path.join(dir_stage, MDFileNames.state)

        with open(state_fn, 'rb') as f:
            state = pickle.load(f)

        return state

    @stage_system
    @mdstages
    def get_stage_topology(self, stg_name='last'):
        """
        This method returns the MD topology of the selected stage name. If no stage name is passed
        the last stage is selected.

        Parameters:
        -----------
        stg_name: String
            The MD stage name

        Return:
        -------
            : OEMol
            The topology of the selected MD stage
        """

        stage = self.get_stage_by_name(stg_name)

        stage_name = stage.get_value(Fields.stage_name)

        dir_stage = self.processed[stage_name]

        topology_fn = os.path.join(dir_stage, MDFileNames.topology)

        topology = oechem.OEMol()

        with oechem.oemolistream(topology_fn) as ifs:
            oechem.OEReadMolecule(ifs, topology)

        return topology

    @stage_system
    @mdstages
    def get_stage_trajectory(self, stg_name='last'):
        """
        This method returns the trajectory file name associated with the md data. If the trajectory is
        not found None is return

        Parameters:
        -----------
        stg_name: String
            The MD stage name

        Return:
        -------
            : String or None
            Trajectory file name if the process was successful otherwise None
        """

        stage = self.get_stage_by_name(stg_name)

        stg_name = stage.get_value(Fields.stage_name)

        stg_type = stage.get_value(Fields.stage_type)

        traj_dir = self.processed[stg_name]

        if not stage.has_field(Fields.trajectory):
            return None

        trj_tar = utils.download_file(stage.get_value(Fields.trajectory), os.path.join(traj_dir, MDFileNames.trajectory))

        with tarfile.open(trj_tar) as tar:
            tar.extractall(path=traj_dir)

        trj_field = stage.get_field(Fields.trajectory.get_name())
        trj_meta = trj_field.get_meta()
        md_engine = trj_meta.get_attribute(Meta.Annotation.Description)

        # TODO I do not like this a lot
        if md_engine == MDEngines.OpenMM and not stg_type == MDStageTypes.FEC:
            traj_fn = glob.glob(os.path.join(traj_dir, '*.h5'))[0]
        elif md_engine == MDEngines.Gromacs:
            traj_fn = glob.glob(os.path.join(traj_dir, '*.xtc'))[0]
        else:  # Yank trajectory
            traj_fn = os.path.join(traj_dir, "experiments/solvent1.nc")

        exists = os.path.isfile(traj_fn)

        if exists:
            return traj_fn
        else:
            raise ValueError("Something went wrong recovering the trajectory")

    def add_new_stage(self,
                      stage_name,
                      stage_type,
                      topology,
                      mdstate,
                      data_fn,
                      append=True,
                      log=None,
                      trajectory_fn=None,
                      trajectory_engine=None,
                      trajectory_orion_ui='OrionFile'):
        """
        This method append a new MD stage to the MD stage record.

        Parameters:
        -----------
        stage_name: String
            The new MD stage name
        stage_type: String
            The MD stage type e.g. SETUP, MINIMIZATION etc.
        topology: OEMol
            The topology
        mdstate: OERecord
            The new mdstate made of state positions, velocities and box vectors
        data_fn: String
            The data file name is used only locally and is linked to the MD data associated
            with the stage. In Orion the data file name is not used
        append: Bool
            If the flag is set to true the stage will be appended to the MD stages otherwise
            the last stage will be overwritten but the new created MD stage
        log: String or None
            Log info
        trajectory: String, Int or None
            The trajectory name or id associated with the new MD stage
        trajectory_engine: String or None
            The MD engine used to generate the new MD stage. Possible names: OpenMM or Gromacs
        trajectory_ui_orion: String
            The trajectory string name to be displayed in the Orion UI
         Return:
        -------
            : Bool
            True if the MD stage creation was successful
        """

        record = OERecord()

        record.set_value(Fields.stage_name, stage_name)
        record.set_value(Fields.stage_type, stage_type)

        if log is not None:
            record.set_value(Fields.log_data, log)

        with TemporaryDirectory() as output_directory:

            top_fn = os.path.join(output_directory, MDFileNames.topology)

            with oechem.oemolostream(top_fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, topology)

            state_fn = os.path.join(output_directory, MDFileNames.state)

            with open(state_fn, 'wb') as f:
                pickle.dump(mdstate, f)

            with tarfile.open(data_fn, mode='w:gz') as archive:
                archive.add(top_fn, arcname=os.path.basename(top_fn))
                archive.add(state_fn, arcname=os.path.basename(state_fn))

        if trajectory_fn is not None:

            if not os.path.isfile(trajectory_fn):
                raise IOError("The trajectory file has not been found: {}".format(trajectory_fn))

            trj_meta = OEFieldMeta()
            trj_meta.set_attribute(Meta.Annotation.Description, trajectory_engine)
            trj_field = OEField(Fields.trajectory.get_name(), Fields.trajectory.get_type(), meta=trj_meta)

        if self.rec.has_field(Fields.md_stages):

            stage_names = self.get_stages_names

            if append:
                if stage_name in stage_names:
                    raise ValueError(
                        "The selected stage name is already present in the MD stages: {}".format(stage_names))

            else:
                if stage_name in stage_names and not stage_name == stage_names[-1]:
                    raise ValueError(
                        "The selected stage name is already present in the MD stages: {}".format(stage_names))

            lf = utils.upload_data(data_fn, collection_id=self.collection_id, shard_name=data_fn)

            record.set_value(Fields.mddata, lf)

            if trajectory_fn is not None:
                lft = utils.upload_file(trajectory_fn, orion_ui_name=trajectory_orion_ui)
                record.set_value(trj_field, lft)

            stages = self.get_stages

            if append:
                stages.append(record)
            else:
                self.delete_stage_by_name('last')
                stages[-1] = record

            self.rec.set_value(Fields.md_stages, stages)

        else:

            lf = utils.upload_data(data_fn, collection_id=self.collection_id, shard_name=data_fn)

            record.set_value(Fields.mddata, lf)

            if trajectory_fn is not None:
                lft = utils.upload_file(trajectory_fn, orion_ui_name=trajectory_orion_ui)
                record.set_value(trj_field, lft)

            self.rec.set_value(Fields.md_stages, [record])

        self.processed[stage_name] = False

        return True

    @property
    @mdstages
    def get_stages(self):
        """
        This method returns the MD stage record list with all the MD stages.

        Return:
        -------
            : list
            The list of the MD stages
        """
        return self.rec.get_value(Fields.md_stages)

    @property
    @mdstages
    def get_stages_names(self):
        """
        This method returns the name list of the MD stages.

        Return:
        -------
            : list
            The list of the MD stage names
        """

        stages = self.rec.get_value(Fields.md_stages)
        stg_names = [stage.get_value(Fields.stage_name) for stage in stages]

        return stg_names

    @property
    def has_stages(self):
        """
        This method returns True if the record has a MD record list, False otherwise.

        Return:
        -------
            : Bool
            True if the record has a list of MD stages, False otherwise
        """
        if not self.rec.has_field(Fields.md_stages):
            return False
        else:
            return True

    @property
    def delete_stages(self):

        stages = self.get_stages

        for stage in stages:
            fid = stage.get_value(Fields.mddata)
            utils.delete_data(fid, collection_id=self.collection_id)

            if stage.get_value(Fields.trajectory) is not None:
                tid = stage.get_value(Fields.trajectory)
                utils.delete_file(tid)

        self.processed = {}
        self.rec.delete_field(Fields.md_stages)

        return True

    def get_parmed(self, sync_stage_name=None):
        """
        This method returns the Parmed object. An exception is raised if the Parmed object cannot
        be found. If sync_stage_name is not None the parmed structure positions, velocities and
        box vectors will be synchronized with the MD State selected by passing the MD stage name

        Parameters:
        -----------
        sync_stage_name: String or None
            The stage name that is used to synchronize the Parmed structure

        Return:
        -------
            : Parmed Structure object
            The Parmed Structure object
        """

        if not self.rec.has_field(Fields.pmd_structure):
            raise ValueError("The Parmed reference is not present on the record")

        pmd_structure = self.rec.get_value(Fields.pmd_structure)

        if in_orion():
            session = APISession

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            collection = session.get_resource(ShardCollection, self.collection_id)

            shard = session.get_resource(Shard(collection=collection), pmd_structure)

            with TemporaryDirectory() as output_directory:

                parmed_fn = os.path.join(output_directory, "parmed.pickle")

                shard.download_to_file(parmed_fn)

                with open(parmed_fn, 'rb') as f:
                    parm_dic = pickle.load(f)

                pmd_structure = parmed.structure.Structure()
                pmd_structure.__setstate__(parm_dic)

        if sync_stage_name is not None:
            mdstate = self.get_stage_state(stg_name=sync_stage_name)

            pmd_structure.positions = mdstate.get_positions()
            pmd_structure.velocities = mdstate.get_velocities()
            pmd_structure.box_vectors = mdstate.get_box_vectors()

        return pmd_structure

    def set_parmed(self, pmdobj, sync_stage_name=None, shard_name=""):
        """
        This method sets the Parmed object. Return True if the setting was successful.
        If sync_stage_name is not None the parmed structure positions, velocities and
        box vectors will be synchronized with the MD State selected by passing the MD
        stage name

        Parameters:
        -----------
        pmjobj: Parmed Structure object
            The Parmed Structure object to be set on the record
        sync_stage_name: String or None
            The stage name that is used to synchronize the Parmed structure
        shard_name: String
            In Orion tha shard will be named by using the shard_name



        Return:
        -------
            : Bool
            True if the setting was successful
        """

        if not isinstance(pmdobj, parmed.Structure):
            raise ValueError("The passed Parmed object is not a valid Parmed Structure: {}".format(type(pmdobj)))

        if sync_stage_name is not None:

            mdstate = self.get_stage_state(stg_name=sync_stage_name)

            pmdobj.positions = mdstate.get_positions()
            pmdobj.velocities = mdstate.get_velocities()
            pmdobj.box_vectors = mdstate.get_box_vectors()

        if in_orion():

            with TemporaryDirectory() as output_directory:

                parmed_fn = os.path.join(output_directory, 'parmed.pickle')

                with open(parmed_fn, 'wb') as f:
                    pickle.dump(pmdobj.__getstate__(), f)

                if self.collection_id is None:
                    raise ValueError("The Collection ID is None")

                if self.rec.has_field(Fields.pmd_structure):
                    fid = self.rec.get_value(Fields.pmd_structure)
                    utils.delete_data(fid, collection_id=self.collection_id)

                session = APISession

                collection = session.get_resource(ShardCollection, self.collection_id)

                shard = Shard.create(collection, name=shard_name)

                shard.upload_file(parmed_fn)

                shard.close()

                self.rec.set_value(Fields.pmd_structure, shard.id)
        else:
            self.rec.set_value(Fields.pmd_structure, pmdobj)

        return True

    @property
    def has_parmed(self):
        """
        This method checks if the Parmed object is on the record.

        Return:
        -------
            : Bool
            True if the Parmed object is on the record otherwise False
        """

        if not self.rec.has_field(Fields.pmd_structure):
            return False
        else:
            return True

    @property
    def delete_parmed(self):
        """
        This method deletes the Parmed object from the record. True is returned if the deletion was
        successful.

        Return:
        -------
            : Bool
            True if Parmed object deletion was successful
        """

        if not self.has_parmed:
            raise ValueError("The Parmed structure is not present on the record")

        if in_orion():

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            session = APISession

            collection = session.get_resource(ShardCollection, self.collection_id)

            file_id = self.rec.get_value(Fields.pmd_structure)

            session.delete_resource(Shard(collection=collection, id=file_id))

        self.rec.delete_field(Fields.pmd_structure)

        return True

    @property
    def get_protein_traj(self):
        """
        This method returns the protein molecule where conformers have been set as trajectory frames


        Return:
        -------
            : OEMol
            The Protein molecule with trajectory conformers
        """

        if not self.rec.has_field(Fields.protein_traj_confs):
            raise ValueError("The protein conformer trajectory is not present on the record")

        protein_conf = self.rec.get_value(Fields.protein_traj_confs)

        if in_orion():

            session = APISession

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            collection = session.get_resource(ShardCollection, self.collection_id)

            shard = session.get_resource(Shard(collection=collection), protein_conf)

            with TemporaryDirectory() as output_directory:

                protein_fn = os.path.join(output_directory, MDFileNames.trajectory_conformers)

                shard.download_to_file(protein_fn)

                protein_conf = oechem.OEMol()

                with oechem.oemolistream(protein_fn) as ifs:
                    oechem.OEReadMolecule(ifs, protein_conf)

        return protein_conf

    def set_protein_traj(self, protein_conf, shard_name=""):
        """
        This method sets the protein trajectory conformers on the record

        Parameters:
        -----------
        protein_conf: OEChem
            The Protein molecule with trajectory conformers
         shard_name: String
            In Orion tha shard will be named by using the shard_name


        Return:
        -------
            : Bool
            True if the setting was successful
        """

        if not isinstance(protein_conf, oechem.OEMol):
            raise ValueError("The passed Parmed object is not a valid Parmed Structure: {}".format(type(protein_conf)))

        if in_orion():

            with TemporaryDirectory() as output_directory:

                protein_fn = os.path.join(output_directory, MDFileNames.trajectory_conformers)

                with oechem.oemolostream(protein_fn) as ofs:
                    oechem.OEWriteConstMolecule(ofs, protein_conf)

                if self.collection_id is None:
                    raise ValueError("The Collection ID is None")

                if self.rec.has_field(Fields.protein_traj_confs):
                    fid = self.rec.get_value(Fields.protein_traj_confs)
                    utils.delete_data(fid, collection_id=self.collection_id)

                session = APISession

                collection = session.get_resource(ShardCollection, self.collection_id)

                shard = Shard.create(collection, name=shard_name)

                shard.upload_file(protein_fn)

                shard.close()

                self.rec.set_value(Fields.protein_traj_confs, shard.id)
        else:
            self.rec.set_value(Fields.protein_traj_confs, protein_conf)

        return True

    def __getattr__(self, name):
        try:
            return getattr(self.rec, name)
        except AttributeError:
            raise AttributeError(
                "'%s' object has no attribute '%s'" % (type(self).__name__, name))
