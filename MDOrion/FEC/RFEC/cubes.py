# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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


from orionplatform.mixins import RecordPortsMixin

from orionplatform.ports import (RecordInputPort,
                                 RecordOutputPort)

from floe.api import ComputeCube

from MDOrion.Standards.standards import Fields

import traceback

from orionplatform.parameters import FileInputParameter

from datarecord.utils import TemporaryPath

from MDOrion.FEC.RFEC import utils

from datarecord import OERecord


class BoundUnboundSwitchCube(RecordPortsMixin, ComputeCube):
    title = "Bound and UnBound Switching Cube"
    # version = "0.1.4"
    classification = [["Simulation Flask Preparation"]]
    tags = ['Simulation', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube emits complexes on the bound port and non-complexes on
    the standard out port. The Flask ids are re-set
    """

    uuid = "80e656e8-33c5-4560-99b5-95e4f69c1701"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    bound_port = RecordOutputPort("bound_port", initializer=False)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count = 0

    def process(self, record, port):

        try:
            if not record.has_value(Fields.md_components):
                raise ValueError("MD Components Field is missing")

            md_components = record.get_value(Fields.md_components)

            record.set_value(Fields.flaskid, self.count)

            if md_components.has_protein:
                record.set_value(Fields.FEC.RBFEC.thd_leg_type, "Bound_OPLMD")
                self.bound_port.emit(record)
            else:
                record.set_value(Fields.FEC.RBFEC.thd_leg_type, "UnBound_OPLMD")
                self.success.emit(record)

            self.count += 1

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class RBFECMapping(RecordPortsMixin, ComputeCube):
    title = "RBFEC Edge Mapping"

    classification = [["Relative Free Energy"]]
    tags = ['Ligand', 'Edge Mapping']
    description = """
    This cube set the edge mapping between the input ligands to perform relative
    binding free energy calculations. The edge mapping is set by providing a text file 
    where each row contains a string with the ligand names involved in the relative 
    binding free energy calculation for example: 

    ....
    lig_i_name >> lig_j_name
    ....

    Input:
    -------
    Data record Stream - Streamed-in of the ligand molecules
    File name - File with the relative binding free energy mapping

    Output:
    -------
    Data Record Stream - Streamed-out of records where ligands have
    been paired to run relative binding free energy calculations.
    """

    uuid = "5b9f7b2f-68e8-4541-a0b7-ddbe6084923f"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    map_file = FileInputParameter("map_file", title="RBFEC Mapping file",
                                  description="RBFEC mapping file", required=True,
                                  default=None)

    bound_port = RecordInputPort("bound_port", initializer=False)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

        self.Bound_edges = dict()
        self.Bound_records = dict()

        self.Unbound_edges = dict()
        self.Unbound_records = dict()

        file = list(self.args.map_file)
        for file_obj in file:
            with TemporaryPath() as path:
                file_obj.copy_to(path)
                with open(path, "r") as f:
                    edge_list = f.readlines()

        if not edge_list:
            raise IOError("Edge file is empty {}")

        for edge in edge_list:

            if not edge:
                continue

            # Comment
            if edge.startswith(";"):
                continue

            edge_str = utils.edge_map_grammar(edge)

            if len(edge_str) > 3:
                raise ValueError("Syntax Error Edge File: {}".format(edge_str))

            edge_A_State = edge_str[0]
            edge_B_State = edge_str[2]

            if edge_A_State == edge_B_State:
                self.opt['Logger'].warn("Edge with the same starting and final state detected. "
                                        "The edge will be skipped")
                continue

            self.count = 0
            self.Bound_edges[(edge_A_State, edge_B_State)] = [False, False]
            self.Unbound_edges[(edge_A_State, edge_B_State)] = [False, False]

    def process(self, record, port):
        try:

            if port == 'intake':
                edges = self.Unbound_edges
                lig_records = self.Unbound_records

            else:
                edges = self.Bound_edges
                lig_records = self.Bound_records

            if not record.has_value(Fields.ligand_name):
                raise ValueError("The record is missing the ligand name field")

            lig_name = record.get_value(Fields.ligand_name)

            # Populate the ligand Bond/Unbond record dictionary
            lig_records[lig_name] = record

            # Mark True all the occurrences of lig_name in the A and B state
            for edge, occurrence in edges.items():

                lig_A_name = edge[0]
                lig_B_name = edge[1]

                if lig_name == lig_A_name:
                    occurrence[0] = True

                elif lig_name == lig_B_name:
                    occurrence[1] = True

            del_edges = [(edge[0], edge[1]) for edge, occ in edges.items()
                         if (occ[0] and occ[1])]

            # print(del_edges)

            for edge_to_del in del_edges:

                if edge_to_del in edges:

                    del edges[edge_to_del]

                    lig_A_name = edge_to_del[0]
                    lig_B_name = edge_to_del[1]

                    lig_A_rec = lig_records[lig_A_name]
                    lig_B_rec = lig_records[lig_B_name]

                    new_record = OERecord()
                    new_record.set_value(Fields.FEC.RBFEC.NESC.state_A, lig_A_rec)
                    new_record.set_value(Fields.FEC.RBFEC.NESC.state_B, lig_B_rec)
                    new_record.set_value(Fields.FEC.RBFEC.edge_name, lig_A_name + '_to_' + lig_B_name)
                    new_record.set_value(Fields.FEC.RBFEC.edgeid, self.count)

                    self.success.emit(new_record)
                    self.count += 1

            lig_rec_to_del = []
            for ln in lig_records.keys():
                to_del = True
                for edge in edges:
                    if ln == edge[0] or ln == edge[1]:
                        to_del = False

                if to_del:
                    lig_rec_to_del.append(ln)

            for ln in lig_rec_to_del:
                del lig_records[ln]

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):

        if self.Bound_edges:
            self.opt['Logger'].warn("The Bound edge list is not empty :{}".format(self.Bound_edges))

        if self.Unbound_edges:
            self.opt['Logger'].warn("The UnBound edge list is not empty :{}".format(self.Unbound_edges))

        return
