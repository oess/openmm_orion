from floe.api.orion import config_from_env
from cuberecord.orion import StreamingDataset
from openeye import oechem
from cuberecord import OERecordComputeCubeBase, OEField
from cuberecord.ports import DataRecordOutputPort
from floe.api import parameter
from floe.constants import ADVANCED
from datarecord import OEReadDataRecord, Types, OEDataRecord, Column
from cuberecord.oldrecordutil import oe_mol_to_data_record
from cuberecord.constants import DEFAULT_MOL_NAME


class ProteinSetReaderCube(OERecordComputeCubeBase):
    success = DataRecordOutputPort('success')
    title = "Record Reader (New Data Model)"
    description = """
    Reads a data set from Orion

    This Reader cube is for the new data model
    """
    classification = [["Input/Output"]]
    tags = ["OpenEye", "OEDataRecord", "Reader", "Input/Output", "Input", "I/O"]

    data_in = parameter.DataSetInputParameter('data_in',
                                              required=True,
                                              title='Dataset to read from',
                                              description='The dataset to read from')

    limit = parameter.IntegerParameter('limit',
                                       required=False,
                                       description='Maximum number of records to read with this cube',
                                       level=ADVANCED)

    chunk_size = parameter.IntegerParameter('chunk_size',
                                            required=False,
                                            default=1000,
                                            min_value=1,
                                            max_value=100000,
                                            description='Number of datarecords to retrieve from the '
                                                        'server for each request',
                                            level=ADVANCED)

    log_timer = parameter.BooleanParameter('log_timer',
                                           title="Enable timing log",
                                           default=False,
                                           description="Log timing of the reader to the log")

    # Enable molecule unique identifier generation
    IDTag = parameter.BooleanParameter('IDTag',
                                       default=True,
                                       required=False,
                                       help_text='If True/Checked proteins are enumerated by sequentially integers.'
                                                 'A data record column is added to the data record with tag: ID')

    # Protein Prefix
    protein_prefix = parameter.StringParameter('protein_prefix',
                                               default='PRT',
                                               help_text='The protein prefix name used to identify the protein')

    def begin(self):
        if self.args.log_timer:
            stopwatch = oechem.OEStopwatch()
        else:
            stopwatch = None

        limit = self.args.limit
        if limit is not None:
            limit = int(limit)

        chunk_size = self.args.chunk_size
        if chunk_size is not None:
            chunk_size = int(chunk_size)

        count = 0
        self.config = config_from_env()
        in_orion = self.config is not None
        if in_orion:
            # Read from Orion's database
            stream = StreamingDataset(self.args.data_in,
                                      input_format='.oeb',
                                      block_size=chunk_size)
            for record in stream:
                # record seems to be an OEDataRecord object not OERecord

                col_mol = Column(DEFAULT_MOL_NAME, Types.Chem.Mol)

                if not col_mol.has_value(record):
                    self.log.warn("Missing '{}' column".format(col_mol.get_name()))
                    self.failure.emit(record)
                    return

                mol = col_mol.get_value(record)
                mol_copy = oechem.OEMol(mol)

                # Try to recognize the residue name
                oechem.OEPerceiveResidues(mol_copy)

                warn = False
                for at in mol_copy.GetAtoms():
                    res = oechem.OEAtomGetResidue(at)
                    if res.GetName() == 'UNL':
                        warn = True
                        break
                if warn:
                    self.log.warn("Unknown residue names detected")

                if self.args.IDTag:
                    name = 'p' + self.args.protein_prefix + '_' + str(count)
                    col_id = Column("ID", Types.String)
                    col_id.set_value(record, name)

                self.success.emit(record)
                count += 1
                if limit is not None and count >= limit:
                    break
        else:
            # Read locally
            if self.args.data_in.endswith(".oedb") or self.args.data_in.endswith(".json"):
                if self.args.data_in.endswith("json"):
                    fmt = 'json'
                else:
                    fmt = 'oeb'
                # Input file is an OEDataRecord binary
                ifs = oechem.oeifstream(str(self.args.data_in))
                while True:
                    try:
                        record = OEReadDataRecord(ifs, fmt=fmt)
                        if record is None:
                            break
                        count += 1
                        if limit is not None and count > limit:
                            break
                        self.success.emit(record)
                    except:
                        ifs.close()
            else:
                # Input is molecules
                with oechem.oemolistream(str(self.args.data_in)) as imstr:
                    for mol in imstr.GetOEMols():

                        mol_copy = oechem.OEMol(mol)
                        record = oe_mol_to_data_record(mol)

                        # Try to recognize the residue name
                        oechem.OEPerceiveResidues(mol_copy)

                        warn = False
                        for at in mol_copy.GetAtoms():
                            res = oechem.OEAtomGetResidue(at)
                            if res.GetName() == 'UNL':
                                warn = True
                                break
                        if warn:
                            self.log.warn("Unknown residue names detected")

                        if self.args.IDTag:
                            name = 'p' + self.args.protein_prefix + '_' + str(count)
                            field_id = OEField("ID", Types.String)
                            record.set_value(field_id, name)

                        if record is None:
                            break
                        count += 1
                        if limit is not None and count > limit:
                            break
                        self.success.emit(record)
        if stopwatch is not None:
            self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
                                                                                 stopwatch.Elapsed(),
                                                                                 count/stopwatch.Elapsed()))