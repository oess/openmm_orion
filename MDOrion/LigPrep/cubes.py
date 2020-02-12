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

from MDOrion.LigPrep import ff_utils

from floe.api import (parameter,
                      ComputeCube,
                      ParallelMixin)

from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)

from oeommtools import utils as oeommutils

from MDOrion.Standards import Fields

from openeye import oechem

from orionplatform.mixins import RecordPortsMixin


class LigandChargeCube(RecordPortsMixin, ComputeCube):
    title = "Ligand Charge"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ["Ligand"]
    description = """
    This cube charges small organic molecules by using the ELF10 charge method 
    (based on am1bcc method). If the ligands are already charged and the user would 
    like to skip this stage the cube parameter “charge_ligand” can be used. 
    The cube requires a record as input with small organic molecules to be charged 
    and produces a new record with the charged molecules.

    Input:
    -------
    oechem.OEMCMol - Streamed-in of molecule to be charged 

    Output:
    -------
    oechem.OEMCMol - Streamed-out of records with the charged molecules.
    """

    uuid = "ea184f6e-feb8-46f1-a89a-6b87270063a3"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers generated to charge the ligands")

    charge_ligands = parameter.BooleanParameter(
        'charge_ligands',
        default=True,
        description='Flag used to set if charge the ligands or not')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing Primary Molecule field")

            ligand = record.get_value(Fields.primary_molecule)

            # Charge the ligand
            if self.opt['charge_ligands']:
                charged_ligand = ff_utils.assignELF10charges(ligand,
                                                             self.opt['max_conformers'],
                                                             strictStereo=False,
                                                             opt=self.opt)

                # If the ligand has been charged then transfer the computed
                # charges to the starting ligand
                map_charges = {at.GetIdx(): at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
                for at in ligand.GetAtoms():
                    at.SetPartialCharge(map_charges[at.GetIdx()])
                self.log.info("[{}] Charges successfully applied to the ligand: {}".format(self.title,
                                                                                           ligand.GetTitle()))

            # Set the primary molecule with the newly charged molecule
            record.set_value(Fields.primary_molecule, ligand)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)


class LigandSetting(RecordPortsMixin, ComputeCube):
    title = "Ligand Setting"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['Ligand']
    description = """
    This cube is used to set the ligand residue name as the cube parameter
    “lig_res_name” (default: “LIG”). This is necessary to facilitate the
    identification of system components during a system splitting.

    Input:
    -------
    Data record Stream - Streamed-in of the ligand molecules

    Output:
    -------
    Data Record Stream - Streamed-out of records where each ligand has
    a new residue name.
    """

    uuid = "fce16dd4-ce3a-4374-92f0-4ed24259d2f6"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    # Ligand Residue Name
    lig_res_name = parameter.StringParameter('lig_res_name',
                                             default='LIG',
                                             help_text='The new ligand residue name')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.ligand_count = 0

    def process(self, initialRecord, port):
        try:
            if not initialRecord.has_value(Fields.primary_molecule):
                raise ValueError("Missing Primary Molecule field")

            ligand = initialRecord.get_value(Fields.primary_molecule)

            # place the entire initial record as a sub-record, to be restored when conformer runs are gathered
            record = OERecord()
            record.set_value(Fields.ligInit_rec, initialRecord)

            if oechem.OECalculateMolecularWeight(ligand) > 900.0:  # Units are in Dalton
                self.opt['Logger'].warn("[{}] The molecule {} seems to have a large molecular "
                                        "weight for a ligand: {:.2f} Da"
                                        .format(self.title,
                                                ligand.GetTitle(),
                                                oechem.OECalculateMolecularWeight(ligand)))

            # Removing Interaction Hint Container, Style and PDB Data
            oechem.OEDeleteInteractionsHintSerializationData(ligand)
            oechem.OEDeleteInteractionsHintSerializationIds(ligand)
            oechem.OEClearStyle(ligand)
            oechem.OEClearPDBData(ligand)

            # Ligand sanitation
            ligand = oeommutils.sanitizeOEMolecule(ligand)

            lig_title = ligand.GetTitle()

            if lig_title == "":
                lig_title = 'LIG'

            record.set_value(Fields.ligand_name, lig_title)

            for at in ligand.GetAtoms():
                residue = oechem.OEAtomGetResidue(at)
                residue.SetName(self.args.lig_res_name)
                oechem.OEAtomSetResidue(at, residue)

            record.set_value(Fields.primary_molecule, ligand)
            record.set_value(Fields.ligid, self.ligand_count)

            self.success.emit(record)
            self.ligand_count += 1

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(initialRecord)


class ParallelLigandChargeCube(ParallelMixin, LigandChargeCube):
    title = "Parallel " + LigandChargeCube.title
    description = "(Parallel) " + LigandChargeCube.description
    uuid = "5f26ca60-301b-4488-9eec-dfb5a760be26"


# class ParallelLigandSetting(ParallelMixin, LigandSetting):
#     title = "Parallel " + LigandSetting.title
#     description = "(Parallel) " + LigandSetting.description
