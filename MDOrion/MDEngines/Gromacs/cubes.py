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

from orionplatform.mixins import RecordPortsMixin

from orionplatform.ports import RecordOutputPort, RecordInputPort

from floe.api import (ParallelMixin,
                      ComputeCube,
                      SourceCube,
                      SinkCube)

from datarecord import OERecord
from datarecord.utils import TemporaryPath

from floe.api import parameters

from orionplatform.parameters import DatasetInputParameter, FileInputParameter

from openeye import oechem

import os

from os import environ

from orionclient.session import in_orion

import tempfile

from orionclient.types import Dataset

from orionclient.session import APISession

from datarecord import OEWriteRecord

from MDOrion.MDEngines.Gromacs.standards import Fields, Gromacs

from MDOrion.MDEngines.Gromacs.utils import gmx_run, gmx_steps


class InputGromacs(SourceCube):
    uuid = "734c4f6f-8ccf-4d78-b37e-8c1a5b143454"
    # version = "0.1.4"
    title = "InputTprGromacs"
    classification = [["Gromacs", "Reader"]]
    tags = ["OpenEye", "Gromacs", "MD"]
    description = "This cube read in a Gromacs .tpr file or a recovery dataset"

    success = RecordOutputPort('success')
    failure = RecordOutputPort('failure')

    prefix_name = parameters.StringParameter("prefix_name",
                                             description="The system prefix name",
                                             required=True,
                                             default="PROT")

    tpr = FileInputParameter("tpr",
                             title="Gromacs Tpr",
                             description="Gromacs 2019 Tpr file input",
                             required=False,
                             default=None)

    data_in = DatasetInputParameter(
        "data_in",
        title="Input Dataset",
        description="The Dataset used for restarting",
        required=False,
        default=None)

    def __iter__(self):
        datasets = list(self.args.data_in)
        files = list(self.args.tpr)
        if len(datasets) > 0:
            for dataset in datasets:
                for record in dataset.records():
                    yield record
        elif len(files) > 0:
            if self.args.prefix_name == "":
                self.args.prefix_name = "PROT"
            for file_obj in files:
                with TemporaryPath(suffix=".tpr") as path:
                    file_obj.copy_to(path)
                    with open(path, "rb") as f:
                        tpr_bytes = f.read()
                        record = OERecord()
                        record.set_value(Fields.tpr_field, tpr_bytes)
                        record.set_value(Fields.prefix_name_field, self.args.prefix_name)
                        yield record
                    break
        else:
            raise ValueError("A Gromacs input .tpr file or a restart dataset is required")


class GromacsProxyCube(RecordPortsMixin, ComputeCube):
    uuid = "d04b093d-f710-4cdc-a62a-f901183c1847"
    # version = "0.1.4"
    title = "Gromacs Proxy Cube"
    description = """
    This cube is used to implement a cycle with the Gromacs running cube
    Checking the current state of the simulation
    """
    classification = [["Proxy"]]
    tags = ["Gromacs", "MD"]

    # Override defaults for some parameters
    parameter_overrides = {
       "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            if record.has_value(Fields.cycle_id):
                current_iteration = record.get_value(Fields.current_iteration_field)
                nsteps = record.get_value(Fields.md_nsteps_field)
            else:
                record.set_value(Fields.cycle_id, 0)
                current_iteration = 0
                record.set_value(Fields.current_iteration_field, 0)

                # Read from the .tpr file the max  number of iterations
                # ad set it on the record
                if not record.has_value(Fields.tpr_field):
                    raise ValueError("Gromacs Tpr field is missing")

                tpr = record.get_value(Fields.tpr_field)

                with tempfile.TemporaryDirectory() as outdir:

                    fn_tpr = os.path.join(outdir, Gromacs.gmx_tpr_prefix + '.tpr')

                    f = open(fn_tpr, 'wb')
                    f.write(tpr)
                    f.close()

                    nsteps = gmx_steps("-s", fn_tpr)

                    record.set_value(Fields.md_nsteps_field, nsteps)

            self.opt['Logger'].info("{} current iterations {}".format(self.title, current_iteration))
            self.opt['Logger'].info("{}  max iterations {}".format(self.title, nsteps))

            if current_iteration == nsteps:
                self.opt['Logger'].info("{} Finishing...".format(self.title))
            else:
                self.opt['Logger'].warn("{} Forwarding to Cycle...".format(self.title))
                self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class GromacsRunCube(RecordPortsMixin, ComputeCube):
    uuid = "50310c8b-4e3a-4f9d-8853-4983363bc247"
    # version = "0.1.4"
    title = "Gromacs Run Cube"
    description = """
    This cube runs Gromacs 
    """
    classification = [["Gromacs"]]
    tags = ["Gromacs", "MD"]

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
       "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    cube_run_time = parameters.DecimalParameter("cube_run_time",
                                                default=10.0,
                                                help_text="Gromacs Max Cube Running Time in hrs")

    verbose = parameters.BooleanParameter("verbose",
                                          default=False,
                                          help_text="Enable/Disable Gromacs Std Output")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        opt = dict(self.opt)

        if not record.has_value(Fields.cycle_id):
            raise ValueError("Missing the current cycle ID")

        cycle_id = record.get_value(Fields.cycle_id)
        opt['cycle_id'] = cycle_id

        if not record.has_value(Fields.tpr_field):
            raise ValueError("Gromacs Tpr field is missing")

        tpr = record.get_value(Fields.tpr_field)
        opt['tpr'] = tpr

        if not record.has_value(Fields.prefix_name_field):
            raise ValueError("System prefix name is missing")

        opt['prefix_name'] = record.get_value(Fields.prefix_name_field)

        opt['record'] = record

        current_md_steps = gmx_run(opt)

        # Update the current iterations and the cycles
        record.set_value(Fields.current_iteration_field, current_md_steps)
        record.set_value(Fields.cycle_id, cycle_id + 1)

        self.success.emit(record)


class WriterRecordCube(SinkCube):
    description = """A cube that takes records and writes them out to a single dataset"""
    title = "Dataset Writer"
    classification = [["I/O", "Writers"]]
    tags = ["I/O", "Writer", "OERecord", "Output"]
    uuid = "7dad346a-70f8-4fd0-b452-7dc0b9b4673c"
    # version = "0.1.4"

    intake = RecordInputPort("intake")

    def write(self, record, port):

        if not record.has_value(Fields.cycle_id):
            raise ValueError("Missing the current cycle ID")

        cycle_id = record.get_value(Fields.cycle_id)

        if not record.has_value(Fields.prefix_name_field):
            raise ValueError("System prefix name is missing")

        prefix_name = record.get_value(Fields.prefix_name_field)

        name_dataset = prefix_name + '_' + "Recovery_Dataset_" + str(cycle_id - 1)

        if in_orion():  # Output to database
            stream = Dataset.create(APISession, name_dataset)
            job_id = environ.get('ORION_JOB_ID')
            APISession.tag_resource(stream, "Job " + str(job_id))
            APISession.tag_resource(stream, "Gmx_Dataset")
            stream.write(record)
            stream.finalize()
        else:
            name_dataset += '.oedb'
            stream = oechem.oeofstream(name_dataset)
            OEWriteRecord(stream, record, fmt="binary")
            stream.close()


# class ParallelGromacsRunCube(ParallelMixin, GromacsRunCube):
#     title = "Parallel " + GromacsRunCube.title
#     description = "(Parallel) " + GromacsRunCube.description
