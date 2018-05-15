from floe.api.parameter import (IntegerParameter,
                                BooleanParameter,
                                StringParameter)

from floe.constants import ADVANCED
from floe.api.orion import in_orion

from orionclient import OrionSession
from orionclient import Dataset

from cuberecord import OERecordSourceCube
from cuberecord.parameters import DataSourceParameter

from cuberecord.cube_testing import OEMolRecordStream
from cuberecord.converters.oldrecordutil import oe_data_record_to_mol

from openeye import oechem

from time import time

from Standards import Fields


# class LigandReaderCube(OERecordSourceCube):
#     title = "Data Set Reader"
#
#     data_in = DataSourceParameter('data_in',
#                                   required=True,
#                                   title='Data to read from',
#                                   description='The data to read from')
#
#     limit = IntegerParameter('limit',
#                              required=False,
#                              description='Maximum number of records to read with this cube',
#                              level=ADVANCED)
#
#     log_timer = BooleanParameter('log_timer',
#                                  title="Enable timing log",
#                                  default=False,
#                                  description="Log timing of the reader to the log")
#
#     # Enable molecule unique identifier generation
#     ID = BooleanParameter('ID',
#                           default=True,
#                           required=False,
#                           help_text='If True/Checked ligands are enumerated by sequentially integers.'
#                                     'An OEField is added to the data record with field: ID',
#                           level=ADVANCED)
#
#     # Ligand Residue Name
#     lig_res_name = StringParameter('lig_res_name',
#                                    default='LIG',
#                                    required=True,
#                                    help_text='The ligand residue name',
#                                    level=ADVANCED)
#
#     def __init__(self, name, **kwargs):
#         super().__init__(name, **kwargs)
#         self._begin_time = None
#         self._count = None
#         self._record_ids = None
#         self._data_set_ids = None
#
#     def begin(self):
#         self._count = 0
#
#         if self.data_in.data_set_key in self.args.data_in:
#             self._data_set_ids = self.args.data_in[self.data_in.data_set_key]
#         else:
#             self._data_set_ids = []
#
#         if self.data_in.record_key in self.args.data_in:
#             self._record_ids = self.args.data_in[self.data_in.record_key]
#         else:
#             self._record_ids = []
#
#         if len(self._record_ids) + len(self._data_set_ids) <= 0:
#             raise ValueError("No input specified")
#
#         if len(self._record_ids) != 0:
#             raise ValueError("Reading record by record ID not supported outside of orion")
#
#         if self.args.log_timer:
#             self._begin_time = time()
#         else:
#             self._begin_time = None
#
#     def __iter__(self):
#         if in_orion():
#             session = OrionSession()
#             for data_set_id in self._data_set_ids:
#                 ifs = session.get_resource(Dataset, data_set_id)
#                 for record in ifs.records():
#                     if self.args.limit is None or self._count < self.args.limit:
#                         self._count += 1
#                         yield record
#                     else:
#                         break
#                 ifs.finalize()
#             if self.data_in.record_key in self.args.data_in:
#                 if len(self.args.data_in[self.data_in.record_key]) == 0:
#                     raise RuntimeError("Reading of individual records by id not supported yet")
#         else:
#             for data_set_id in self._data_set_ids:
#                 ifs = OEMolRecordStream(data_set_id)
#                 for record in ifs:
#                     #mol = record.get_value(Fields.primary_molecule)
#
#                     mol = oe_data_record_to_mol(record)
#
#                     for at in mol.GetAtoms():
#                         residue = oechem.OEAtomGetResidue(at)
#                         residue.SetName(self.args.lig_res_name)
#                         oechem.OEAtomSetResidue(at, residue)
#
#                     if self.args.ID:
#                         name = 'l' + mol.GetTitle()[0:12] + '_' + str(self._count)
#                         record.set_value(Fields.id, name)
#
#                     record.set_value(Fields.primary_molecule, mol)
#
#                     if self.args.limit is None or self._count < self.args.limit:
#                         self._count += 1
#                         yield record
#                     else:
#                         break
#
#     def end(self):
#         if self._begin_time is not None:
#             run_time = time() - self._begin_time
#             self.log.info("Read {0} records in {1} second.  {2}records/second".format(self._count,
#                                                                                       run_time,
#                                                                                       self._count/run_time))
#

# class LigandReaderCube(SourceCube):
#     success = RecordOutputPort('success')
#     title = "Ligand Reader"
#     description = """
#     Reads a data set from Orion
#
#     This Ligand Reader cube is for the new data model
#     """
#     classification = [["Input/Output"]]
#     tags = ["OpenEye", "OERecord", "Reader", "Input/Output", "Input", "I/O"]
#
#     data_in = StringParameter('data_in',
#                               required=True,
#                               title='Data set to read from',
#                               description='The data set to read from')
#
#     limit = IntegerParameter('limit',
#                              required=False,
#                              description='Maximum number of records to read with this cube',
#                              level=ADVANCED)
#
#     chunk_size = IntegerParameter('chunk_size',
#                                   required=False,
#                                   default=1000,
#                                   min_value=1,
#                                   max_value=100000,
#                                   description='Number of Records to retrieve from the server for each request',
#                                   level=ADVANCED)
#     log_timer = BooleanParameter('log_timer',
#                                  title="Enable timing log",
#                                  default=False,
#                                  description="Log timing of the reader to the log")
#
#     # Enable molecule unique identifier generation
#     ID = BooleanParameter('ID',
#                           default=True,
#                           required=False,
#                           help_text='If True/Checked ligands are enumerated by sequentially integers.'
#                                     'An OEField is added to the data record with field: ID',
#                           level=ADVANCED)
#
#     # Ligand Residue Name
#     lig_res_name = StringParameter('lig_res_name',
#                                    default='LIG',
#                                    required=True,
#                                    help_text='The ligand residue name',
#                                    level=ADVANCED)
#
#     def begin(self):
#         if self.args.log_timer:
#             self._stopwatch = oechem.OEStopwatch()
#             self._begin_time = time()
#         else:
#             self._stopwatch = None
#             self._begin_time = None
#         self.limit = self.args.limit
#         if self.limit is not None:
#             self.limit = int(self.limit)
#
#         self.chunk_size = self.args.chunk_size
#         if self.chunk_size is not None:
#             self.chunk_size = int(self.chunk_size)
#
#         self.config = config_from_env()
#         in_orion = self.config is not None
#
#         self.stream = None
#         self.ifs = None
#         self.imstr = None
#
#         if in_orion:
#             # Read from Orion's database
#             self.stream = StreamingDataset(self.args.data_in,
#                                            input_format='.oeb',
#                                            block_size=self.chunk_size)
#         else:
#             # Read locally
#             if self.args.data_in.endswith(".oedb") or self.args.data_in.endswith(".json"):
#                 if self.args.data_in.endswith("json"):
#                     self.fmt = 'json'
#                 else:
#                     self.fmt = 'oeb'
#                 # Input file is an OERecord binary
#                 self.ifs = oechem.oeifstream(str(self.args.data_in))
#             else:
#                 # Input is molecules
#                 self.imstr = oechem.oemolistream(str(self.args.data_in))
#
#     def __iter__(self):
#         count = 0
#         if self.stream:
#             for record in self.stream:
#                 yield record
#                 count += 1
#                 if self.limit is not None and self.count >= self.limit:
#                     break
#         elif self.ifs:
#             while True:
#                 try:
#                     record = OEReadDataRecord(self.ifs, fmt=self.fmt)
#                     if record is None:
#                         break
#                     count += 1
#                     if self.limit is not None and self.count > self.limit:
#                         break
#                     yield record
#                 # except:
#                 #     self.ifs.close()
#                 finally:
#                     self.ifs.close()
#         elif self.imstr:
#             for mol in self.imstr.GetOEMols():
#
#                 for at in mol.GetAtoms():
#                     residue = oechem.OEAtomGetResidue(at)
#                     residue.SetName(self.args.lig_res_name)
#                     oechem.OEAtomSetResidue(at, residue)
#
#                 record = oe_mol_to_data_record(mol, include_sd_data=False)
#
#                 if self.args.ID:
#                     name = 'l' + mol.GetTitle()[0:12] + '_' + str(count)
#                     record.set_value(Fields.id, name)
#
#                 if record is None:
#                     break
#                 count += 1
#                 if self.limit is not None and self.count > self.limit:
#                     break
#                 yield record
#             self.imstr.close()
#         else:
#             self.log.error('No input method found')
#
#         if self._stopwatch is not None:
#             elapsed = self._stopwatch.Elapsed()
#             self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
#                                                                                  elapsed,
#                                                                                  count/elapsed))
#         if self._begin_time is not None:
#             wall_time = time() - self._begin_time
#             self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
#                                                                                  wall_time,
#                                                                                  count/wall_time))


# class LigandReaderCube(OERecordComputeCubeBase):
#     success = DataRecordOutputPort('success')
#     title = "Record Reader (New Data Model)"
#     description = """
#     Reads a data set from Orion
#
#     This Reader cube is for the new data model
#     """
#     classification = [["Input/Output"]]
#     tags = ["OpenEye", "OEDataRecord", "Reader", "Input/Output", "Input", "I/O"]
#
#     data_in = parameter.DataSetInputParameter('data_in',
#                                               required=True,
#                                               title='Dataset to read from',
#                                               description='The dataset to read from')
#
#     limit = parameter.IntegerParameter('limit',
#                                        required=False,
#                                        description='Maximum number of records to read with this cube',
#                                        level=ADVANCED)
#
#     chunk_size = parameter.IntegerParameter('chunk_size',
#                                             required=False,
#                                             default=1000,
#                                             min_value=1,
#                                             max_value=100000,
#                                             description='Number of datarecords to retrieve from the server '
#                                                         'for each request',
#                                             level=ADVANCED)
#
#     log_timer = parameter.BooleanParameter('log_timer',
#                                            title="Enable timing log", default=False,
#                                            description="Log timing of the reader to the log")
#
#     # Enable molecule unique identifier generation
#     IDTag = parameter.BooleanParameter('IDTag',
#                                        default=True,
#                                        required=False,
#                                        help_text='If True/Checked ligands are enumerated by sequentially integers.'
#                                                  'A data record column is added to the data record with tag: IDTag')
#
#     # Molecule Type
#     type_lig = parameter.StringParameter('type_lig',
#                                          default='LIG',
#                                          required=True,
#                                          help_text='The ligand residue name. A data record column is added to'
#                                                    ' the data record with tag: type_lig')
#
#     def begin(self):
#         if self.args.log_timer:
#             stopwatch = oechem.OEStopwatch()
#         else:
#             stopwatch = None
#
#         limit = self.args.limit
#         if limit is not None:
#             limit = int(limit)
#
#         chunk_size = self.args.chunk_size
#         if chunk_size is not None:
#             chunk_size = int(chunk_size)
#
#         count = 0
#         self.config = config_from_env()
#         in_orion = self.config is not None
#         if in_orion:
#             # Read from Orion's database
#             stream = StreamingDataset(self.args.data_in,
#                                       input_format='.oeb',
#                                       block_size=chunk_size)
#             for record in stream:
#
#                 col_ligand = Column(DEFAULT_MOL_NAME, Types.Chem.Mol)
#
#                 if not col_ligand.has_value(record):
#                     self.log.warn("Missing '{}' column".format(col_ligand.get_name()))
#                     self.failure.emit(record)
#                     return
#
#                 ligand = col_ligand.get_value(record)
#
#                 for at in ligand.GetAtoms():
#                     residue = oechem.OEAtomGetResidue(at)
#                     residue.SetName(self.args.type_lig)
#                     oechem.OEAtomSetResidue(at, residue)
#
#                 col_ligand.set_value(record, ligand)
#
#                 if self.args.IDTag:
#                     name = 'l' + ligand.GetTitle()[0:12] + '_' + str(count)
#                     col_id = Column("ID", Types.String)
#                     col_id.set_value(record, name)
#
#                 self.success.emit(record)
#                 count += 1
#                 if limit is not None and count >= limit:
#                     break
#         else:
#             # Read locally
#             if self.args.data_in.endswith(".oedb"):
#                 # Input file is an OEDataRecord binary
#                 ifs = oechem.oeifstream(str(self.args.data_in))
#                 while True:
#                     try:
#                         record = OEReadDataRecord(ifs, fmt='oeb')
#                         if record is None:
#                             break
#                         count += 1
#                         if limit is not None and count > limit:
#                             break
#                         self.success.emit(record)
#                     except:
#                         ifs.close()
#             else:
#                 # Input is molecules
#                 with oechem.oemolistream(str(self.args.data_in)) as imstr:
#                     for mol in imstr.GetOEMols():
#
#                         for at in mol.GetAtoms():
#                             residue = oechem.OEAtomGetResidue(at)
#                             residue.SetName(self.args.type_lig)
#                             oechem.OEAtomSetResidue(at, residue)
#
#                         record = oe_mol_to_data_record(mol, include_sd_data=False)
#
#                         if self.args.IDTag:
#                             name = 'l' + mol.GetTitle()[0:12] + '_' + str(count)
#                             field_id = OEField("ID", Types.String)
#                             record.set_value(field_id, name)
#
#                         if record is None:
#                             break
#                         count += 1
#                         if limit is not None and count > limit:
#                             break
#                         self.success.emit(record)
#         if stopwatch is not None:
#             self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
#                                                                                  stopwatch.Elapsed(),
#                                                                                  count/stopwatch.Elapsed()))
#
#
# class DataSetWriterCubeStripCustom(OERecordComputeCubeBase):
#     intake = DataRecordInputPort('intake')
#     title = "Record Writer (New Data Model)"
#     description = """
#     Writes a data set to Orion
#
#     This writer cube is for the new data model"""
#     classification = [["Input/Output"]]
#     tags = ["OpenEye", "OEDataRecord", "Reader", "Input/Output", "Output", "I/O"]
#
#     data_out = parameter.DataSetOutputParameter('data_out',
#                                                 required=True,
#                                                 title='Name of data set to create',
#                                                 description='The data set to output')
#
#     chunk_size = parameter.IntegerParameter('chunk_size',
#                                             required=False,
#                                             default=1000,
#                                             min_value=1,
#                                             max_value=10000,
#                                             description='Number of data records to send to the server with each request',
#                                             level=ADVANCED)
#
#     log_timer = parameter.BooleanParameter('log_timer',
#                                            title="Enable timing log",
#                                            default=False,
#                                            description="Log timing of the reader to the log")
#
#     pack_mol = parameter.BooleanParameter('pack_mol',
#                                           default=True,
#                                           title='Pack Molecule',
#                                           description='Pack Child Molecular data for the orion backend.  '
#                                                       'Only turn this off if '
#                                                       'you are doing this explicity in an earlier cube.')
#
#     _stopwatch = None
#     _count = 0
#     _begin_time = None
#     _process_time = None
#     config = None
#     in_orion = None
#     ofs = None
#
#     def begin(self):
#
#         if self.args.log_timer:
#             self._stopwatch = oechem.OEStopwatch()
#             self._count = 0
#             self._begin_time = time()
#             self._process_time = 0.0
#         self.config = config_from_env()
#         self.in_orion = self.config is not None
#         if self.in_orion:
#             if self.args.log_timer:
#                 log = self.log
#             else:
#                 log = None
#             self.ofs = SimpleDatasetUploader(self.args.data_out,
#                                              chunk_size=self.args.chunk_size,
#                                              tags=[self.name],
#                                              log=log)
#         else:
#             if self.args.data_out.endswith(".oedb"):
#                 self.ofs = oechem.oeofstream(str(self.args.data_out))
#             else:
#                 self.ofs = oechem.oemolostream(str(self.args.data_out))
#
#     def process(self, record, port):
#         if self._process_time is not None:
#             begin_time = time()
#         else:
#             begin_time = None
#         if self.in_orion:
#             # col_Parmed = Column("Parmed", Types.Custom)
#             #
#             # if col_Parmed.has_value(record):
#             #     col_Parmed.delete_column(record)
#
#             for col in record.get_columns(include_meta=False):
#                 if col.get_type() == Types.Custom:
#                     col.delete_from(record)
#
#             self.ofs.write_record(record)
#         else:
#
#             # for col in record.get_columns(include_meta=False):
#             #     if col.get_type() == Types.Custom:
#             #         print(col)
#             #         col.delete_from(record)
#
#             if self.args.data_out.endswith(".oedb"):
#                 OEWriteDataRecord(self.ofs, record, fmt='oeb')
#             elif self.args.data_out.endswith(".json"):
#                 OEWriteDataRecord(self.ofs, record, fmt='json')
#             else:
#                 mol = oe_data_record_to_mol(record)
#                 if mol is not None:
#                     oechem.OEWriteMolecule(self.ofs, mol)
#         self._count += 1
#         if self._process_time is not None:
#             self._process_time += time() - begin_time
#
#     def end(self):
#         self.ofs.close()
#         if self._stopwatch is not None:
#             self.log.info("Wrote {} molecules in {} seconds. (CPU time)"
#                           "({} molecules/sec".format(self._count,
#                                                      self._stopwatch.Elapsed(),
#                                                      self._count/self._stopwatch.Elapsed()))
#         if self._process_time is not None and self._process_time > 0.0:
#             self.log.info("Process function Wrote {} molecules in {} seconds (wall time). "
#                           "({} molecules/sec".format(self._count,
#                                                      self._process_time,
#                                                      self._count/self._process_time))
#         if self._begin_time is not None:
#             wall_time = time() - self._begin_time
#             self.log.info("Wrote {} molecules in {} process seconds. (wall time)"
#                           "({} molecules/sec".format(self._count,
#                                                      wall_time,
#                                                      self._count/wall_time))
