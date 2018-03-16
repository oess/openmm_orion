import traceback
from LigPrepCubes import ff_utils
from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube, OEField
from datarecord import Types, Meta, ColumnMeta
from oeommtools import utils as oeommutils
from cuberecord.constants import DEFAULT_MOL_NAME

from tempfile import TemporaryDirectory
import os
import random
import string
import tarfile
from LigPrepCubes.ff_utils import upload, download
from floe.api.orion import in_orion, upload_file
from big_storage import LargeFileDataType
from OpenMMCubes.utils import ParmedData


class LigandSetChargeCube(ParallelMixin, OERecordComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
    This cube charges the Ligand by using the ELF10 charge method

    Input:
    -------
    oechem.OEMCMol - Streamed-in of the ligand molecules

    Output:
    -------
    oechem.OEMCMol - Emits the charged ligands
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers")

    charge_ligands = parameter.BooleanParameter(
        'charge_ligands',
        default=True,
        description='Flag used to set if charge the ligands or not')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            field_mol = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)

            if not record.has_value(field_mol):
                self.log.warn("Missing '{}' field".format(field_mol.get_name()))
                self.failure.emit(record)
                return

            ligand = record.get_value(field_mol)

            # Ligand sanitation
            ligand = oeommutils.sanitizeOEMolecule(ligand)

            # Charge the ligand
            if self.opt['charge_ligands']:
                charged_ligand = ff_utils.assignELF10charges(ligand,
                                                             self.opt['max_conformers'],
                                                             strictStereo=False)

                # If the ligand has been charged then transfer the computed
                # charges to the starting ligand
                map_charges = {at.GetIdx(): at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
                for at in ligand.GetAtoms():
                    at.SetPartialCharge(map_charges[at.GetIdx()])
                self.log.info("ELF10 charge method applied to the ligand: {}".format(ligand.GetTitle()))

            record.set_value(field_mol, ligand)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            self.log.warn("Failed to assign ELF10 charges on molecule {}".format(ligand.GetTitle()))
            # Return failed record
            self.failure.emit(record)