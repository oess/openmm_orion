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
                                 OrionSession,
                                 get_session)

import glob


from orionclient.helpers.collections import (try_hard_to_create_shard,
                                             try_hard_to_download_shard)

from oemdtoolbox.ForceField.md_components import MDComponents


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
    """

    This Class Implements the MD Datarecord API by using
    getter and setter functions

    """

    def __init__(self, record, inplace=True):
        """
        The Initialization function used to create the
        MDDatarecord object

        Parameters
        ---------
        record: OERecord object
            The OERecord used to create the MDDatarecord
        inplace: Bool
            if True the record will be update in place otherwise
            a copy of the record will be made
        """
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

        Parameters
        ---------

        Returns
        -------
        record: OERecord
            The record to be passed with the cubes
        """
        return self.rec

    @property
    def get_primary(self):
        """
        This method returns the primary molecule present on the record

        Parameters
        ---------

        Returns
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

        Parameters
        -----------

        primary_mol: OEMol
            The primary molecule to set on the record

        Returns
        -------
        boolean: Bool
            True if the primary molecule has been set on the record
        """

        if not isinstance(primary_mol, oechem.OEMol):
            raise ValueError("The Primary Molecule is not a valid OEMol: {}".format(primary_mol))

        self.rec.set_value(Fields.primary_molecule, primary_mol)

        return True

    @property
    def get_ligand(self):
        """
        This method returns the ligand molecule present on the record

        Parameters
        -----------

        Returns
        -------
        ligand: OEMol
            The ligand molecule if the ligand has been set on the record otherwise
            an error is raised
        """

        if not self.rec.has_field(Fields.ligand):
            raise ValueError("The ligand molecule has not been found on the record")

        return self.rec.get_value(Fields.ligand)

    @property
    def has_ligand(self):
        """
        This method returns True if ligand molecule is present on the record

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
            True if the ligand molecule is present on the record otherwise False
        """

        if self.rec.has_field(Fields.ligand):
            return True
        else:
            return False

    def set_ligand(self, ligand):
        """
        This method sets the ligand molecule  on the record

        Parameters
        ----------

        Returns
        -------
        boolean: Bool
            returns True if the ligand has been set on the record otherwise an error
            is raised
        """

        if not isinstance(ligand, oechem.OEMol):
            raise ValueError("The ligand molecule is not a valid OEMol: {}".format(ligand))

        self.rec.set_value(Fields.ligand, ligand)

        return True

    @property
    def get_protein(self):
        """
        This method returns the protein molecule present on the record

        Parameters
        ----------

        Returns
        -------
        protein : OEMol
            The protein molecule if the protein has been set on the record otherwise
            an error is raised
        """

        if not self.rec.has_field(Fields.protein):
            raise ValueError("The protein molecule has not been found on the record")

        return self.rec.get_value(Fields.protein)

    @property
    def has_protein(self):
        """
        This method returns true if the protein molecule is present on the record

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
            True if the protein molecule is present on the record otherwise False
        """

        if self.rec.has_field(Fields.protein):
            return True
        else:
            return False

    def set_protein(self, protein):
        """
        This method sets the protein molecule  on the record

        Parameters
        ----------

        Returns
        -------
        boolean: Bool
            returns True if the protein has been set on the record otherwise an
            error is raised
        """

        if not isinstance(protein, oechem.OEMol):
            raise ValueError("The protein molecule is not a valid OEMol: {}".format(protein))

        self.rec.set_value(Fields.protein, protein)

        return True

    @property
    def get_flask(self):
        """
        This method returns the flask molecule present on the record

        Parameters
        ----------

        Returns
        -------

        flask: OEMol
            The flask present on the record otherwise an error is raised
        """

        if not self.rec.has_field(Fields.flask):
            raise ValueError("The Flask Molecule has not been found on the record")

        return self.rec.get_value(Fields.flask)

    def set_flask(self, flask):
        """
        This method sets the flask molecule on the record

        Parameters
        ----------
        flask : OEMol
            The flask molecule to set on the record

        Returns
        -------
        record: Bool
            True if the flask molecule has been set on the record
            otherwise an error is raised
        """

        if not isinstance(flask, oechem.OEMol):
            raise ValueError("The flask Molecule is not a valid OEMol: {}".format(flask))

        self.rec.set_value(Fields.flask, flask)

        return True

    @property
    def get_flask_id(self):
        """
        This method returns the integer value of the flask identification field present on the record

        Parameters
        ----------

        Returns
        -------
        flask_id : Int
            The unique flask identifier
        """

        if not self.rec.has_field(Fields.flaskid):
            raise ValueError("The flask ID Field has not been found on the record")

        return self.rec.get_value(Fields.flaskid)

    @property
    def has_flask_id(self):
        """
        This method checks if the flask identification field is present on the record

        Parameters
        ----------

        Returns
        -------
        boolean: Bool
            True if the flask ID field is present on the record otherwise False
        """

        if not self.rec.has_field(Fields.flaskid):
            return False
        else:
            return True

    def set_flask_id(self, id):
        """
        This method sets the integer value of the flask identification field on the record

        Parameters
        -----------
        id: Int
            An integer value for the flask identification field

        Returns
        -------
        boolean : Bool
            True if the flask identification ID has been set as an integer on the record
        """

        if not isinstance(id, int):
            raise ValueError(" The id must be an integer: {}".format(id))

        self.rec.set_value(Fields.flaskid, id)

        return True

    @property
    def get_lig_id(self):
        """
        This method returns the ligand identification field present on the record

        Parameters
        -----------

        Returns
        -------
        ligand_id: Int
            The ligand identification id number
        """

        if not self.rec.has_field(Fields.ligid):
            raise ValueError("The ligand identification field has not been found on the record")

        return self.rec.get_value(Fields.ligid)

    @property
    def has_lig_id(self):
        """
        This method checks if the ligand identification field is present on the record

        Parameters
        -----------

        Returns
        -------
        boolean : Bool
            True if the ligand identification field is present on the record otherwise False
        """

        if not self.rec.has_field(Fields.ligid):
            return False
        else:
            return True

    def set_lig_id(self, sys_id):
        """
        This method sets the ligand identification field on the record

        Parameters
        ----------
        sys_id: Int
            An integer value for the ligand identification field

        Returns
        -------
        boolean: Bool
            True if the value for the ligand identification field was successfully set on the record
        """

        if not isinstance(sys_id, int):
            raise ValueError("The system id must be an integer: {}".format(sys_id))

        self.rec.set_value(Fields.ligid, sys_id)

        return True

    @property
    def get_conf_id(self):
        """
        This method returns the identification field CONF ID present on the record

        Parameters
        ----------

        Returns
        -------
        conformed_id: Int
            The conformer id
        """

        if not self.rec.has_field(Fields.confid):
            raise ValueError("The CONF ID Field has not been found on the record")

        return self.rec.get_value(Fields.confid)

    @property
    def has_conf_id(self):
        """
        This method checks if the identification field CONF ID is present on the record

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
            True if the conformer id field is present on the record otherwise False
        """

        if not self.rec.has_field(Fields.confid):
            return False
        else:
            return True

    def set_conf_id(self, conf_id):
        """
        This method sets the identification field for the conformer on the record

        Parameters
        -----------
        conf_id: Int
            An identification integer for the record

        Returns
        -------
        boolean : Bool
            True if the conformed id has been set on the record otherwise an error is raised
        """

        if not isinstance(conf_id, int):
            raise ValueError("The conformer id must be an integer: {}".format(conf_id))

        self.rec.set_value(Fields.confid, conf_id)

        return True

    @property
    def get_title(self):
        """
        This method returns the title present on the record

        Parameters
        -----------

        Returns
        -------
        title : String
            The title string if present on the record otherwise an error is raised
        """

        if not self.rec.has_field(Fields.title):
            raise ValueError("The Title Field has not been found on the record")

        return self.rec.get_value(Fields.title)

    @property
    def has_title(self):
        """
        This method checks if the Title field is present on the record

        Parameters
        -----------

        Returns
        -------
        boolean : Bool
            True if the Title field is resent on the record otherwise False
        """

        if not self.rec.has_field(Fields.title):
            return False
        else:
            return True

    def set_title(self, title):
        """
        This method sets the system Title field on the record

        Parameters
        ----------
        title: String
            A string used to identify the molecular system

        Returns
        -------
        boolean : Bool
            True if the system Title has been set on the record
        """

        if not isinstance(title, str):
            raise ValueError(" The title must be a sting: {}".format(title))

        self.rec.set_value(Fields.title, title)

        return True

    def create_collection(self, name):
        """
        This method sets a collection field on the record to be used in Orion

        Parameters
        -----------
        name: String
            A string used to identify in the Orion UI the collection

        Returns
        -------
        boolean : Bool
                True if the collection creation in Orion was successful otherwise False
        """

        if in_orion():

            if self.rec.has_field(Fields.collection):
                raise ValueError("Collection field already present on the record")

            # session = APISession

            session = OrionSession(
                requests_session=get_session(
                    retry_dict={
                        403: 5,
                        404: 20,
                        409: 45,
                        460: 15,
                        500: 2,
                        502: 45,
                        503: 45,
                        504: 45,
                    }
                )
            )

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

        Parameters
        -----------

        Returns
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

        Parameters
        -----------
        stg_name: String
            The MD stage name


        Returns
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
        This method returns a MD stage selected by passing an index. If the stage is not found
        an exception is raised.

        Parameters
        -----------
        idx: Int
            The stage index to retrieve

        Returns
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

        Parameters
        -----------
        stg_name: String
            The MD stage name

        Returns
        -------
        boolean: Bool
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

        Parameters
        -----------
        idx: Int
            The MD stage index

        Returns
        -------
        boolean: Bool
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

        Parameters
        -----------
        stg_name: String
            The MD stage name

        Returns
        -------
        boolean: Bool
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

        Parameters
        ----------
        stg_name: String
            The MD stage name

        Returns
        -------
        info_string: String
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

        Parameters
        -----------
        stg_name: String
            The MD stage name

        Returns
        -------
        state : MDState
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

        Parameters
        -----------
        stg_name: String
            The MD stage name

        Returns
        -------
        topology : OEMol
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
        This method returns the trajectory file name associated
        with the md data. If the trajectory is not found None is return

        Parameters
        ----------
        stg_name: String
            The MD stage name

        Returns
        -------
        trajectory_filename: String or None
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

        if md_engine == MDEngines.OpenMM and not stg_type == MDStageTypes.FEC:
            traj_fn = glob.glob(os.path.join(traj_dir, '*.h5'))[0]
        elif md_engine == MDEngines.Gromacs:
            traj_fn = glob.glob(os.path.join(traj_dir, '*.trr'))[0]
        else:
            raise ValueError("MD Engine Not Supported")

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
        This method add a new MD stage to the MD stage record

        Parameters
        ----------
        stage_name: String
            The new MD stage name
        stage_type: String
            The MD stage type e.g. SETUP, MINIMIZATION etc.
        topology: OEMol
            The topology
        mdstate: MDState
            The new mdstate made of state positions, velocities and box vectors
        data_fn: String
            The data file name is used only locally and is linked to the MD data associated
            with the stage. In Orion the data file name is not used
        append: Bool
            If the flag is set to true the stage will be appended to the MD stages otherwise
            the last stage will be overwritten by the new created MD stage
        log: String or None
            Log info
        trajectory_fn: String, Int or None
            The trajectory name for local run or id in Orion associated with the new MD stage
        trajectory_engine: String or None
            The MD engine used to generate the new MD stage. Possible names: OpenMM or Gromacs
        trajectory_orion_ui: String
            The trajectory string name to be displayed in the Orion UI

        Returns
        -------
        boolean: Bool
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

        Parameters
        ----------

        Returns
        -------
        record_list: list
            The MD stages record list
        """
        return self.rec.get_value(Fields.md_stages)

    @property
    @mdstages
    def get_stages_names(self):
        """
        This method returns the list names of the MD stages.

        Parameters
        ----------

        Returns
        -------
        list: list
            The MD stage name list
        """

        stages = self.rec.get_value(Fields.md_stages)
        stg_names = [stage.get_value(Fields.stage_name) for stage in stages]

        return stg_names

    @property
    def has_stages(self):
        """
        This method returns True if the record has a MD record list otherwise False

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
            True if the record has a list of MD stages otherwise False
        """
        if not self.rec.has_field(Fields.md_stages):
            return False
        else:
            return True

    @property
    def delete_stages(self):
        """
        This method deletes all the record stages

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
            True if the deletion was successful
        """

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

        Parameters
        ----------
        sync_stage_name: String or None
            The stage name that is used to synchronize the Parmed structure

        Returns
        -------
        parmed : Parmed Structure
            The Parmed Structure object
        """

        if not self.rec.has_field(Fields.pmd_structure):
            raise ValueError("The Parmed reference is not present on the record")

        pmd_structure = self.rec.get_value(Fields.pmd_structure)

        if in_orion():
            # session = APISession

            session = OrionSession(
                requests_session=get_session(
                    retry_dict={
                        403: 5,
                        404: 20,
                        409: 45,
                        460: 15,
                        500: 2,
                        502: 45,
                        503: 45,
                        504: 45,
                    }
                )
            )

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            collection = session.get_resource(ShardCollection, self.collection_id)

            shard = session.get_resource(Shard(collection=collection), pmd_structure)

            with TemporaryDirectory() as output_directory:

                parmed_fn = os.path.join(output_directory, "parmed.pickle")

                try_hard_to_download_shard(shard, parmed_fn)

                with open(parmed_fn, 'rb') as f:
                    parm_dic = pickle.load(f)

                pmd_structure = parmed.structure.Structure()
                pmd_structure.__setstate__(parm_dic)

            shard.close()

        if sync_stage_name is not None:
            mdstate = self.get_stage_state(stg_name=sync_stage_name)

            pmd_structure.positions = mdstate.get_positions()
            pmd_structure.velocities = mdstate.get_velocities()
            pmd_structure.box_vectors = mdstate.get_box_vectors()

        return pmd_structure

    def set_parmed(self, pmd, sync_stage_name=None, shard_name=""):
        """
        This method sets the Parmed object. Return True if the setting was successful.
        If sync_stage_name is not None the parmed structure positions, velocities and
        box vectors will be synchronized with the MD State selected by passing the MD
        stage name

        Parameters
        ----------
        pmd: Parmed Structure object
            The Parmed Structure object to be set on the record
        sync_stage_name: String or None
            The stage name that is used to synchronize the Parmed structure
        shard_name: String
            In Orion tha shard will be named by using the shard_name

        Returns
        -------
        boolean : Bool
            True if the setting was successful
        """

        if not isinstance(pmd, parmed.Structure):
            raise ValueError("The passed Parmed object is not a valid Parmed Structure: {}".format(type(pmd)))

        if sync_stage_name is not None:

            mdstate = self.get_stage_state(stg_name=sync_stage_name)

            pmd.positions = mdstate.get_positions()
            pmd.velocities = mdstate.get_velocities()
            pmd.box_vectors = mdstate.get_box_vectors()

        if in_orion():

            with TemporaryDirectory() as output_directory:

                parmed_fn = os.path.join(output_directory, 'parmed.pickle')

                with open(parmed_fn, 'wb') as f:
                    pickle.dump(pmd.__getstate__(), f)

                if self.collection_id is None:
                    raise ValueError("The Collection ID is None")

                if self.rec.has_field(Fields.pmd_structure):
                    fid = self.rec.get_value(Fields.pmd_structure)
                    utils.delete_data(fid, collection_id=self.collection_id)

                session = OrionSession(
                    requests_session=get_session(
                        retry_dict={
                            403: 5,
                            404: 20,
                            409: 45,
                            460: 15,
                            500: 2,
                            502: 45,
                            503: 45,
                            504: 45,
                        }
                    )
                )

                collection = session.get_resource(ShardCollection, self.collection_id)

                shard = try_hard_to_create_shard(collection, parmed_fn, name=shard_name)

                shard.close()

                self.rec.set_value(Fields.pmd_structure, shard.id)
        else:
            self.rec.set_value(Fields.pmd_structure, pmd)

        return True

    @property
    def has_parmed(self):
        """
        This method checks if the Parmed object is on the record.

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
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

        Parameters
        ----------

        Returns
        -------
        boolean : Bool
            True if Parmed object deletion was successful
        """

        if not self.has_parmed:
            raise ValueError("The Parmed structure is not present on the record")

        if in_orion():

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            session = OrionSession(
                requests_session=get_session(
                    retry_dict={
                        403: 5,
                        404: 20,
                        409: 45,
                        460: 15,
                        500: 2,
                        502: 45,
                        503: 45,
                        504: 45,
                    }
                )
            )

            collection = session.get_resource(ShardCollection, self.collection_id)

            file_id = self.rec.get_value(Fields.pmd_structure)

            session.delete_resource(Shard(collection=collection, id=file_id))

        self.rec.delete_field(Fields.pmd_structure)

        return True

    @property
    def get_protein_traj(self):
        """
        This method returns the protein molecule where conformers have been set as trajectory frames

        Parameters
        ----------
        Returns
        -------
        multi_conformer_protein: OEMol
            The multi conformer protein
        """

        if not self.rec.has_field(Fields.protein_traj_confs):
            raise ValueError("The protein conformer trajectory is not present on the record")

        protein_conf = self.rec.get_value(Fields.protein_traj_confs)

        if in_orion():

            # session = APISession

            session = OrionSession(
                requests_session=get_session(
                    retry_dict={
                        403: 5,
                        404: 20,
                        409: 45,
                        460: 15,
                        500: 2,
                        502: 45,
                        503: 45,
                        504: 45,
                    }
                )
            )

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            collection = session.get_resource(ShardCollection, self.collection_id)

            shard = session.get_resource(Shard(collection=collection), protein_conf)

            with TemporaryDirectory() as output_directory:

                protein_fn = os.path.join(output_directory, MDFileNames.trajectory_conformers)

                try_hard_to_download_shard(shard, protein_fn)

                protein_conf = oechem.OEMol()

                with oechem.oemolistream(protein_fn) as ifs:
                    oechem.OEReadMolecule(ifs, protein_conf)

            shard.close()

        return protein_conf

    def set_protein_traj(self, protein_conf, shard_name=""):
        """
        This method sets the multi conformer protein trajectory on the record

        Parameters
        -----------
        protein_conf: OEChem
            Th multi conformer protein trajectory
        shard_name: String
            In Orion tha shard will be named by using the shard_name

        Returns
        -------
        boolean: Bool
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

                # session = APISession

                session = OrionSession(
                    requests_session=get_session(
                        retry_dict={
                            403: 5,
                            404: 20,
                            409: 45,
                            460: 15,
                            500: 2,
                            502: 45,
                            503: 45,
                            504: 45,
                        }
                    )
                )

                collection = session.get_resource(ShardCollection, self.collection_id)

                shard = try_hard_to_create_shard(collection, protein_fn, name=shard_name)

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

    @property
    def get_md_components(self):
        """
        This method returns the MD Components if present on the record

        Parameters
        -----------

        Returns
        -------
        md_components: MDComponents
            The MD Components object
        """

        if not self.rec.has_field(Fields.md_components):
            raise ValueError("The MD Components Field has not been found on the record")

        return self.rec.get_value(Fields.md_components)

    def set_md_components(self, md_components):
        """
        This method sets the MD Components on the record

        Parameters
        ----------
        md_components: MDComponents
            The MD Components instance

        Returns
        -------
        boolean : Bool
            True if the the md components field was successfully set on the record
        """

        if not isinstance(md_components, MDComponents):
            raise ValueError("The passed components is not an object of type MDComponents: {}".format(type(md_components)))

        self.rec.set_value(Fields.md_components, md_components)

        return True

    @property
    def has_md_components(self):
        """
        This method returns True if the MD Components object is present on the record

        Return:
        -------
        boolean: Bool
            True if the md components object is present on the record otherwise False
        """

        if self.rec.has_field(Fields.md_components):
            return True
        else:
            return False

    @property
    def has_extra_data_tar(self):
        """
        This method returns True if extra data file in tar format is attached to the record

        Return:
        -------
        boolean: Bool
            True if extra data file is present on the record otherwise False
        """

        if self.rec.has_field(Fields.extra_data_tar):
            return True
        else:
            return False

    @property
    def get_extra_data_tar(self):
        """
        This method returns the directory file name where extra data file tar has been unpacked

        Parameters
        ----------
        Returns
        -------
        directory file name: String
            The directory file name
        """

        if not self.rec.has_field(Fields.extra_data_tar):
            raise ValueError("Extra data file tar has not been found on the record")

        extra_data = self.rec.get_value(Fields.extra_data_tar)

        if in_orion():

            # session = APISession

            session = OrionSession(
                requests_session=get_session(
                    retry_dict={
                        403: 5,
                        404: 20,
                        409: 45,
                        460: 15,
                        500: 2,
                        502: 45,
                        503: 45,
                        504: 45,
                    }
                )
            )

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            collection = session.get_resource(ShardCollection, self.collection_id)

            shard = session.get_resource(Shard(collection=collection), extra_data)

            extra_data = os.path.join(self.cwd, 'extra_data.tar')

            try_hard_to_download_shard(shard, extra_data)

            shard.close()

        return extra_data

    def set_extra_data_tar(self, tar_fn, shard_name=""):
        """
        This method sets the extra data file on the record

        Parameters
        -----------
        tar_fn: String
            The compressed data file name
        shard_name: String
            In Orion tha shard will be named by using the shard_name

        Returns
        -------
        boolean: Bool
            True if the setting was successful
        """

        if not os.path.isfile(tar_fn):
            raise ValueError("The filename does not exist: {}".format(tar_fn))

        if in_orion():

            if self.collection_id is None:
                raise ValueError("The Collection ID is None")

            if self.rec.has_field(Fields.extra_data_tar):
                fid = self.rec.get_value(Fields.extra_data_tar)
                utils.delete_data(fid, collection_id=self.collection_id)

            # session = APISession
            session = OrionSession(
                requests_session=get_session(
                    retry_dict={
                        403: 5,
                        404: 20,
                        409: 45,
                        460: 15,
                        500: 2,
                        502: 45,
                        503: 45,
                        504: 45,
                    }
                )
            )

            collection = session.get_resource(ShardCollection, self.collection_id)

            shard = try_hard_to_create_shard(collection, tar_fn, name=shard_name)

            shard.close()

            self.rec.set_value(Fields.extra_data_tar, shard.id)
        else:
            self.rec.set_value(Fields.extra_data_tar, tar_fn)

        return True
