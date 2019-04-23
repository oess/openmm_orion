# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


from datarecord.types import CustomHandler

import pickle

import parmed

from MDOrion.MDEngines.utils import MDState

import copy

from orionclient.session import in_orion, APISession

from orionclient.types import File

from os import environ

import os

from orionclient.types import (Shard,
                               ShardCollection)


class ParmedData(CustomHandler):

    @staticmethod
    def get_name():
        return 'Parmed'

    @classmethod
    def validate(cls, value):
        return isinstance(value, parmed.structure.Structure)

    @classmethod
    def copy(cls, value):
        return parmed.structure.copy(value)

    @staticmethod
    def serialize(structure):
        struct_dict = structure.__getstate__()
        pkl_obj = pickle.dumps(struct_dict)
        return bytes(pkl_obj)

    @staticmethod
    def deserialize(data):
        new_structure = parmed.structure.Structure()
        new_structure.__setstate__(pickle.loads(bytes(data)))
        return new_structure


class MDStateData(CustomHandler):

    @staticmethod
    def get_name():
        return 'MDState'

    @classmethod
    def validate(cls, value):
        return isinstance(value, MDState)

    @classmethod
    def copy(cls, value):
        return copy.deepcopy(value)

    @staticmethod
    def serialize(state):
        pkl_obj = pickle.dumps(state)
        return bytes(pkl_obj)

    @staticmethod
    def deserialize(data):
        new_state = pickle.loads(bytes(data))
        return new_state


def upload_file(filename, orion_ui_name='OrionFile'):

    if in_orion():

        session = APISession

        file_upload = File.upload(session,
                                  orion_ui_name,
                                  filename)

        session.tag_resource(file_upload, "Trajectory")

        job_id = environ.get('ORION_JOB_ID')

        if job_id:
            session.tag_resource(file_upload, "Job {}".format(job_id))

        file_id = file_upload.id

    else:
        file_id = filename

    return file_id


def download_file(file_id, filename):

    if in_orion():

        session = APISession

        resource = session.get_resource(File, file_id)

        resource.download_to_file(filename)

        fn_local = filename

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The trajectory file has not been found: {}".format(fn_local))

    return fn_local


def delete_file(file_id):

    if in_orion():

        session = APISession

        resource = session.get_resource(File, file_id)

        session.delete_resource(resource)

    else:
        os.remove(file_id)

    return True


def upload_data(filename, collection_id=None, shard_name=""):

    if in_orion():

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        session = APISession

        collection = session.get_resource(ShardCollection, collection_id)

        shard = Shard.create(collection, name=shard_name)

        shard.upload_file(filename)

        shard.close()

        file_id = shard.id

    else:
        file_id = filename

    return file_id


def download_data(file_id, path, collection_id=None):

    if in_orion():

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        session = APISession

        collection = session.get_resource(ShardCollection, collection_id)

        shard = session.get_resource(Shard(collection=collection), file_id)

        from MDOrion.Standards import MDFileNames

        fn_local = os.path.join(path, MDFileNames.mddata)

        shard.download_to_file(fn_local)

    else:
        fn_local = file_id

    if not os.path.isfile(fn_local):
        raise IOError("The MD data file has not been found: {}".format(fn_local))

    return fn_local


def delete_data(file_id, collection_id=None):

    if in_orion():

        if collection_id is None:
            raise ValueError("The Collection ID is None")

        session = APISession

        collection = session.get_resource(ShardCollection, collection_id)

        session.delete_resource(Shard(collection=collection, id=file_id))

    else:
        os.remove(file_id)

    return True
