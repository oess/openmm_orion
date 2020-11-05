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

from floe.api import (parameters,
                      ComputeCube)

from MDOrion.Standards.standards import Fields

import traceback

from orionplatform.parameters import FileInputParameter

from datarecord.utils import TemporaryPath

from MDOrion.FEC.RFEC import utils

from datarecord import OERecord

from oemdtoolbox.FEC.RBFEC.chimera import Chimera

from MDOrion.FEC.RFEC.utils import gmx_chimera_topology_injection

from MDOrion.FEC.RFEC.utils import parmed_find_ligand

from MDOrion.Standards.mdrecord import MDDataRecord


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
    # version = "0.1.4"
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

    uuid = "e079546e-a9e7-4d1b-875e-cd180b785992"

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

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.A_State = list()
        self.B_State = list()
        self.ligand_names = list()
        self.ligand_dic = dict()

        file = list(self.args.map_file)
        for file_obj in file:
            with TemporaryPath() as path:
                file_obj.copy_to(path)
                with open(path, "r") as f:
                    map_list = f.readlines()

        if not map_list:
            raise IOError("Map file is empty {}")

        for m in map_list:

            if not m:
                continue

            # Comment
            if m.startswith(";"):
                continue

            expr = utils.edge_map_grammar(m)

            if len(expr) > 3:
                raise ValueError("Syntax Error Map File: {}".format(expr))

            self.A_State.append(expr[0])
            self.B_State.append(expr[2])

    def process(self, record, port):
        try:
            lig_name = record.get_value(Fields.ligand_name)

            if lig_name in self.ligand_names:
                raise ValueError("All ligands must have different names. Duplicate: {}".format(lig_name))
            else:
                self.ligand_names.append(lig_name)

            self.ligand_dic[lig_name] = record

            return

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

    def end(self):
        try:
            state_A_set = set(self.A_State)
            state_B_set = set(self.B_State)
            state_AB_set = state_A_set.union(state_B_set)

            if not state_AB_set:
                raise ValueError("The provide map will not produce any edge with the provided ligands")

            for lig_name in state_AB_set:
                if lig_name in self.ligand_dic:
                    rec = self.ligand_dic[lig_name]
                    self.success.emit(rec)
                else:
                    raise ValueError("The following ligand name has not been found: {}".format(lig_name))

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())


class RBFECEdgeGathering(RecordPortsMixin, ComputeCube):
    title = "RBFEC Edge Gathering"

    classification = [["Relative Free Energy"]]
    tags = ['Ligand', 'Edge Mapping']
    description = """
    TBD
    """

    uuid = "970052ea-25dc-4aa6-87ee-840506022a8b"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
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

        self.Bound_Unbound_edges = dict()

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
            self.Bound_Unbound_edges[(edge_A_State, edge_B_State)] = [self.Bound_edges, self.Unbound_edges]

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

            # Populate the ligand Bound/Unbound record dictionary
            lig_records[lig_name] = record

            # Mark True all the occurrences of lig_name in the A and B state
            for edge, occurrence in edges.items():

                lig_A_name = edge[0]
                lig_B_name = edge[1]

                if lig_name == lig_A_name:
                    occurrence[0] = True

                elif lig_name == lig_B_name:
                    occurrence[1] = True

            full_edges_to_del = [(edge[0], edge[1]) for edge, v in self.Bound_Unbound_edges.items() if v[0][edge][0]
                                 and v[0][edge][1] and v[1][edge][0] and v[1][edge][1]]

            for edge_to_del in full_edges_to_del:
                if edge_to_del in edges:

                    del self.Bound_edges[edge_to_del]
                    del self.Unbound_edges[edge_to_del]

                    lig_A_name = edge_to_del[0]
                    lig_B_name = edge_to_del[1]

                    bound_rec_A = self.Bound_records[lig_A_name]
                    bound_rec_B = self.Bound_records[lig_B_name]

                    unbound_rec_A = self.Unbound_records[lig_A_name]
                    unbound_rec_B = self.Unbound_records[lig_B_name]

                    new_record = OERecord()
                    new_record.set_value(Fields.FEC.RBFEC.NESC.state_A, [bound_rec_A, unbound_rec_A])
                    new_record.set_value(Fields.FEC.RBFEC.NESC.state_B, [bound_rec_B, unbound_rec_B])
                    new_record.set_value(Fields.FEC.RBFEC.edge_name, lig_A_name + '_to_' + lig_B_name)
                    new_record.set_value(Fields.FEC.RBFEC.edgeid, self.count)
                    self.success.emit(new_record)
                    self.count += 1

            # Clean Bound and Unbound Records
            bound_rec_to_del = list()
            for ln in self.Bound_records.keys():
                to_del = True
                for ed in self.Bound_edges:
                    if ln == ed[0] or ln == ed[1]:
                        to_del = False
                if to_del:
                    bound_rec_to_del.append(ln)

            for ln in bound_rec_to_del:
                del self.Bound_records[ln]

            unbound_rec_to_del = list()
            for ln in self.Unbound_records.keys():
                to_del = True
                for ed in self.Unbound_edges:
                    if ln == ed[0] or ln == ed[1]:
                        to_del = False
                if to_del:
                    unbound_rec_to_del.append(ln)

            for ln in unbound_rec_to_del:
                del self.Unbound_records[ln]

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


class GMXChimera(RecordPortsMixin, ComputeCube):
    title = "GMX Chimera"

    classification = [["Relative Free Energy"]]
    tags = ['Ligand', 'Edge Mapping']
    description = """
    TBD
    """

    uuid = "0ccfad69-bada-48fa-a890-348d7784d7ea"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 32000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    # Ligand Residue Name
    lig_res_name = parameters.StringParameter('lig_res_name',
                                              default='LIG',
                                              help_text='The new ligand residue name')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            rec_list_state_A = record.get_value(Fields.FEC.RBFEC.NESC.state_A)
            rec_list_state_B = record.get_value(Fields.FEC.RBFEC.NESC.state_B)

            print(record.get_value(Fields.FEC.RBFEC.edge_name))

            md_record_state_A_Bound = MDDataRecord(rec_list_state_A[0])
            md_record_state_B_Bound = MDDataRecord(rec_list_state_B[0])

            md_record_state_A_Unbound = MDDataRecord(rec_list_state_A[1])
            md_record_state_B_Unbound = MDDataRecord(rec_list_state_B[1])

            lig_A = md_record_state_A_Unbound.get_md_components.get_ligand
            lig_B = md_record_state_B_Unbound.get_md_components.get_ligand

            pmd_flask_state_A_Unbound = md_record_state_A_Unbound.get_parmed(sync_stage_name="System Parametrization")
            pmd_flask_state_B_Unbound = md_record_state_B_Unbound.get_parmed(sync_stage_name="System Parametrization")

            pmd_flask_state_A_Bound = md_record_state_A_Bound.get_parmed(sync_stage_name="System Parametrization")
            pmd_flask_state_B_Bound = md_record_state_B_Bound.get_parmed(sync_stage_name="System Parametrization")

            pmd_lig_A, idxA = parmed_find_ligand(pmd_flask_state_A_Unbound)
            pmd_lig_B, idxB = parmed_find_ligand(pmd_flask_state_B_Unbound)

            if pmd_lig_A is None:
                raise ValueError("It was not possible to extract the ligand parmed from the state A")

            if pmd_lig_B is None:
                raise ValueError("It was not possible to extract the ligand parmed from the state B")

            chimera = Chimera(lig_A, lig_B, pmd_lig_A, pmd_lig_B)

            pmd_chimera_A_to_B_initial, pmd_chimera_A_to_B_final = chimera.pmd_chimera(morph="A_to_B")
            pmd_chimera_B_to_A_initial, pmd_chimera_B_to_A_final = chimera.pmd_chimera(morph="B_to_A")

            gmx_A_to_B_Unbound = gmx_chimera_topology_injection(pmd_flask_state_A_Unbound,
                                                                pmd_chimera_A_to_B_initial,
                                                                pmd_chimera_A_to_B_final)

            gmx_B_to_A_Unbound = gmx_chimera_topology_injection(pmd_flask_state_B_Unbound,
                                                                pmd_chimera_B_to_A_initial,
                                                                pmd_chimera_B_to_A_final)

            gmx_A_to_B_Bound = gmx_chimera_topology_injection(pmd_flask_state_A_Bound,
                                                              pmd_chimera_A_to_B_initial,
                                                              pmd_chimera_A_to_B_final)

            gmx_B_to_A_Bound = gmx_chimera_topology_injection(pmd_flask_state_B_Bound,
                                                              pmd_chimera_B_to_A_initial,
                                                              pmd_chimera_B_to_A_final)

            import sys
            sys.exit(-1)



            self.success.emit(rec_list_state_A[0])

            del md_record_state_A_Bound
            del md_record_state_B_Bound
            del md_record_state_A_Unbound
            del md_record_state_B_Unbound

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)








# class RBFECEdgeGathering(RecordPortsMixin, ComputeCube):
#     title = "RBFEC Edge Gathering"
#
#     classification = [["Relative Free Energy"]]
#     tags = ['Ligand', 'Edge Mapping']
#     description = """
#     TBD
#     """
#
#     uuid = "5b9f7b2f-68e8-4541-a0b7-ddbe6084923f"
#
#     # Override defaults for some parameters
#     parameter_overrides = {
#         "memory_mb": {"default": 14000},
#         "spot_policy": {"default": "Prohibited"},
#         "prefetch_count": {"default": 1},  # 1 molecule at a time
#         "item_count": {"default": 1}  # 1 molecule at a time
#     }
#
#     map_file = FileInputParameter("map_file", title="RBFEC Mapping file",
#                                   description="RBFEC mapping file", required=True,
#                                   default=None)
#
#     bound_port = RecordInputPort("bound_port", initializer=False)
#
#     def begin(self):
#         self.opt = vars(self.args)
#         self.opt['Logger'] = self.log
#
#         self.Bound_edges = dict()
#         self.Bound_records = dict()
#
#         self.Unbound_edges = dict()
#         self.Unbound_records = dict()
#
#         file = list(self.args.map_file)
#         for file_obj in file:
#             with TemporaryPath() as path:
#                 file_obj.copy_to(path)
#                 with open(path, "r") as f:
#                     edge_list = f.readlines()
#
#         if not edge_list:
#             raise IOError("Edge file is empty {}")
#
#         for edge in edge_list:
#
#             if not edge:
#                 continue
#
#             # Comment
#             if edge.startswith(";"):
#                 continue
#
#             edge_str = utils.edge_map_grammar(edge)
#
#             if len(edge_str) > 3:
#                 raise ValueError("Syntax Error Edge File: {}".format(edge_str))
#
#             edge_A_State = edge_str[0]
#             edge_B_State = edge_str[2]
#
#             if edge_A_State == edge_B_State:
#                 self.opt['Logger'].warn("Edge with the same starting and final state detected. "
#                                         "The edge will be skipped")
#                 continue
#
#             self.count = 0
#             self.Bound_edges[(edge_A_State, edge_B_State)] = [False, False]
#             self.Unbound_edges[(edge_A_State, edge_B_State)] = [False, False]
#
#     def process(self, record, port):
#         try:
#
#             if port == 'intake':
#                 edges = self.Unbound_edges
#                 lig_records = self.Unbound_records
#
#             else:
#                 edges = self.Bound_edges
#                 lig_records = self.Bound_records
#
#             if not record.has_value(Fields.ligand_name):
#                 raise ValueError("The record is missing the ligand name field")
#
#             lig_name = record.get_value(Fields.ligand_name)
#
#             # Populate the ligand Bond/Unbond record dictionary
#             lig_records[lig_name] = record
#
#             # Mark True all the occurrences of lig_name in the A and B state
#             for edge, occurrence in edges.items():
#
#                 lig_A_name = edge[0]
#                 lig_B_name = edge[1]
#
#                 if lig_name == lig_A_name:
#                     occurrence[0] = True
#
#                 elif lig_name == lig_B_name:
#                     occurrence[1] = True
#
#             del_edges = [(edge[0], edge[1]) for edge, occ in edges.items()
#                          if (occ[0] and occ[1])]
#
#             # print(del_edges)
#
#             for edge_to_del in del_edges:
#
#                 if edge_to_del in edges:
#
#                     del edges[edge_to_del]
#
#                     lig_A_name = edge_to_del[0]
#                     lig_B_name = edge_to_del[1]
#
#                     lig_A_rec = lig_records[lig_A_name]
#                     lig_B_rec = lig_records[lig_B_name]
#
#                     new_record = OERecord()
#                     new_record.set_value(Fields.FEC.RBFEC.NESC.state_A, lig_A_rec)
#                     new_record.set_value(Fields.FEC.RBFEC.NESC.state_B, lig_B_rec)
#                     new_record.set_value(Fields.FEC.RBFEC.edge_name, lig_A_name + '_to_' + lig_B_name)
#                     new_record.set_value(Fields.FEC.RBFEC.edgeid, self.count)
#
#                     self.success.emit(new_record)
#                     self.count += 1
#
#             lig_rec_to_del = []
#             for ln in lig_records.keys():
#                 to_del = True
#                 for edge in edges:
#                     if ln == edge[0] or ln == edge[1]:
#                         to_del = False
#
#                 if to_del:
#                     lig_rec_to_del.append(ln)
#
#             for ln in lig_rec_to_del:
#                 del lig_records[ln]
#
#         except Exception as e:
#
#             print("Failed to complete", str(e), flush=True)
#             self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
#             self.log.error(traceback.format_exc())
#             self.failure.emit(record)
#
#         return
#
#     def end(self):
#
#         if self.Bound_edges:
#             self.opt['Logger'].warn("The Bound edge list is not empty :{}".format(self.Bound_edges))
#
#         if self.Unbound_edges:
#             self.opt['Logger'].warn("The UnBound edge list is not empty :{}".format(self.Unbound_edges))
#
#         return
