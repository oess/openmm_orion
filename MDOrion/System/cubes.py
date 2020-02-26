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

import traceback

from MDOrion.Standards import Fields

from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.System.utils import get_human_readable

from openeye import oechem

from orionplatform.mixins import RecordPortsMixin

from orionplatform.ports import (RecordInputPort,
                                 RecordOutputPort)


from floe.api import (ParallelMixin,
                      parameter,
                      ComputeCube)

from oeommtools import packmol


from orionclient.types import ShardCollection

from orionclient.session import (in_orion,
                                 APISession)

from os import environ


class IDSettingCube(RecordPortsMixin, ComputeCube):
    title = "Simulation Flask ID Setting"
    # version = "0.1.4"
    classification = [["Simulation Flask Preparation"]]
    tags = ['Simulation', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube sets the integer ID for each simulation flask as well as a descriptive
    title string. If the input molecule 
    on a record has multiple conformers these are split into singles each with 
    its own ID. If a complex will be formed, this cube should be used on ligands
    before forming the complex.
    
    Input:
    -------
    Data record Stream - Streamed input of ligands, one per record

    Output:
    -------
    Data record Stream - Streamed output of records, one per conformer, with title and ID.
    """

    uuid = "d3c1dac4-544f-4273-8b17-1b75c058f4bd"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.total_count = 0
        self.ligid = -1

    def process(self, record, port):
        try:
            if not record.has_value(Fields.flask):
                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("Primary Molecule is missing")
                flask = record.get_value(Fields.primary_molecule)

                record.set_value(Fields.flask, flask)

            flask = record.get_value(Fields.flask)

            # There should be a ligid; if not, increment the last one
            if not record.has_value(Fields.ligid):
                self.ligid += 1
                record.set_value(Fields.ligid, self.ligid)

            if flask.NumConfs() > 1:
                self.opt['Logger'].info("[{}] The flask {} has multiple conformers. Each single conformer "
                                        "will be treated as a new molecule".format(self.title,
                                                                                   flask.GetTitle()))

            name = flask.GetTitle()[0:12]
            if not name:
                name = 'SYS'

            num_conf_counter = 0
            for conf in flask.GetConfs():

                conf_mol = oechem.OEMol(conf)

                flask_title = name

                if flask.GetMaxConfIdx() > 1:
                    flask_title += '_c' + str(num_conf_counter)

                conf_mol.SetTitle(flask_title)

                record.set_value(Fields.flaskid, self.total_count)
                record.set_value(Fields.confid, num_conf_counter)
                record.set_value(Fields.title, flask_title)
                record.set_value(Fields.flask, conf_mol)

                num_conf_counter += 1

                self.total_count += 1
                self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class CollectionSetting(RecordPortsMixin, ComputeCube):
    title = "Collection Setting"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['System', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube set a record collection state in open or closed for safety by
    using the cube bool parameter open. A True value will open the record
    collection enabling the shard writing and deleting. If on the record
    the collection field is not present one will be created.

    Input:
    -------
    Data record Stream - Streamed-in of systems such as ligands

    Output:
    -------
    Data Record Stream - Streamed-out of records each one with associated IDs
    """

    uuid = "b3821952-a5ed-4028-867c-3f71185442aa"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    open = parameter.BooleanParameter(
        'open',
        default=True,
        help_text='Open or Close a Collection')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.collection = None

    def process(self, record, port):
        try:

            if in_orion():

                session = APISession

                if record.has_value(Fields.collection):

                    if self.collection is None:

                        collection_id = record.get_value(Fields.collection)

                        collection = session.get_resource(ShardCollection, collection_id)

                        self.collection = collection

                        if self.opt['open']:

                            self.collection.open()

                else:
                    if self.collection is None:

                        job_id = environ.get('ORION_JOB_ID')

                        self.collection = ShardCollection.create(session, job_id)

                        job_id = environ.get('ORION_JOB_ID')

                        if job_id:
                            session.tag_resource(self.collection, "Job {}".format(job_id))

                    record.set_value(Fields.collection, self.collection.id)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):
        if in_orion():
            if not self.opt['open']:
                if self.collection is not None:
                    self.collection.close()


class SolvationCube(RecordPortsMixin, ComputeCube):
    title = "Solvation Packmol"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['Complex', 'Protein', 'Ligand', 'Solvation']
    description = """
    The solvation cube solvates a given solute input system in a
    selected mixture of solvents. The solvents can be specified by
    comma separated smiles strings of each solvent component or
    selected keywords like tip3p for tip3p water geometry. For each
    component the user needs to specify its molar fractions as well.
    The solution can be neutralized by adding counter-ions. In addition,
    the ionic solution strength can be set adding salt. The cube
    requires a record as input with a solute molecule to solvate
    and produces an output record with the solvated solute.


     Input:
    -------
    Data record Stream - Streamed-in of system solutes to solvate

    Output:
    -------
    Data Record Stream - Streamed-out of records each with the solvated
    solute
    """

    uuid = "2e6130f6-2cba-48a4-9ef3-351a2970258a"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    density = parameter.DecimalParameter(
        'density',
        default=1.0,
        help_text="Solution density in g/ml")

    padding_distance = parameter.DecimalParameter(
        'padding_distance',
        default=8.0,
        help_text="The padding distance between the solute and the box edge in A")

    distance_between_atoms = parameter.DecimalParameter(
        'distance_between_atoms',
        default=2.0,
        help_text="The minimum distance between atoms in A")

    solvents = parameter.StringParameter(
        'solvents',
        default='tip3p',
        help_text='Select solvents. The solvents are specified as comma separated smiles strings'
                  'e.g. [H]O[H], C(Cl)(Cl)Cl, CS(=O)C or special keywords like tip3p')

    molar_fractions = parameter.StringParameter(
        'molar_fractions',
        default='1.0',
        help_text="Molar fractions of each solvent components. The molar fractions are specified"
                  "as comma separated molar fractions strings e.g. 0.5,0.2,0.3")

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text='Output Packmol log')

    geometry = parameter.StringParameter(
        'geometry',
        default='box',
        choices=['box', 'sphere'],
        help_text="Geometry selection: box or sphere. Sphere cannot be used as periodic system "
                  "along with MD simulation")

    close_solvent = parameter.BooleanParameter(
        'close_solvent',
        default=False,
        help_text="If Checked/True solvent molecules will be placed very close to the solute")

    salt = parameter.StringParameter(
        'salt',
        default='[Na+], [Cl-]',
        help_text='Salt type. The salt is specified as list of smiles strings. '
                  'Each smiles string is the salt component dissociated in the '
                  'solution e.g. Na+, Cl-')

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=0.0,
        help_text="Salt concentration in millimolar")

    neutralize_solute = parameter.BooleanParameter(
        'neutralize_solute',
        default=True,
        help_text='Neutralize the solute by adding Na+ and Cl- counter-ions based on'
                  'the solute formal charge')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = dict(self.opt)

            if not record.has_value(Fields.flask):
                raise ValueError("Missing the Flask Molecule Field")

            solute = record.get_value(Fields.flask)

            if not record.has_value(Fields.title):
                self.log.warn("Missing Title field")
                solute_title = solute.GetTitle()[0:12]
            else:
                solute_title = record.get_value(Fields.title)

            self.log.info("[{}] solvating flask {}".format(self.title, solute_title))

            # Update cube simulation parameters
            for field in record.get_fields(include_meta=True):
                field_name = field.get_name()
                if field_name in ['molar_fractions', 'density', 'solvents']:
                    rec_value = record.get_value(field)
                    if field_name == 'molar_fractions':
                        opt[field_name] = str(rec_value)
                    else:
                        opt[field_name] = rec_value
                    opt['Logger'].info("{} Updating parameters for molecule: {} {} = {}".format(self.title,
                                                                                                solute.GetTitle(),
                                                                                                field_name,
                                                                                                rec_value))
            # Solvate the system
            sol_system = packmol.oesolvate(solute, **opt)

            self.log.info("[{}] Solvated simulation flask {} yielding {} atoms overall".format(self.title,
                                                                                              solute_title,
                                                                                              sol_system.NumAtoms()))
            sol_system.SetTitle(solute.GetTitle())

            record.set_value(Fields.flask, sol_system)
            record.set_value(Fields.title, solute_title)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class RecordSizeCheck(RecordPortsMixin, ComputeCube):
    title = "Record Size Checking"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['System', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube checks if the size of the incoming record is less than 100MB
    to avoid Orion database size issues. Locally does not have any effect.

    Input:
    -------
    Data record Stream - Streamed-in of system records

    Output:
    -------
    Data Record Stream - Streamed-out of records
    """

    uuid = "0555ead8-0339-41f2-9876-3eb166e32772"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    fail_in = RecordInputPort("fail_in", initializer=False)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            if in_orion():
                # Create the MD record to use the MD Record API
                mdrecord = MDDataRecord(record)

                if not mdrecord.has_title:
                    self.log.warn("Missing record Title field")
                    system = mdrecord.get_flask
                    title = system.GetTitle()[0:12]
                else:
                    title = mdrecord.get_title

                tot_size = 0
                for field in record.get_fields():
                    tot_size += record.get_value_size(field)

                if tot_size > 100 * 1024 * 1024:
                    raise ValueError("The record size exceeds the 100 MB: {} = {}".format(title,
                                                                                          get_human_readable(tot_size)))
                else:
                    self.opt['Logger'].info("Record size: {} = {}".format(title, get_human_readable(tot_size)))

            if port == "intake":
                self.success.emit(record)
            else:  # Fail in port
                self.failure.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())

        return


class ParallelSolvationCube(ParallelMixin, SolvationCube):
    title = "Parallel " + SolvationCube.title
    description = "(Parallel) " + SolvationCube.description
    uuid = "568ffd29-23e0-4d35-b37c-727596bedf92"


class ParallelRecordSizeCheck(ParallelMixin, RecordSizeCheck):
    title = "Parallel " + RecordSizeCheck.title
    description = "(Parallel) " + RecordSizeCheck.description
    uuid = "f93acfba-a9e8-482b-bcf7-e181e6cb6b09"

