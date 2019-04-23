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
from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube

from oeommtools import utils as oeommutils

from MDOrion.Standards import Fields

from openeye import oechem


class LigandChargeCube(ParallelMixin, OERecordComputeCube):
    title = "Ligand Charge"
    version = "0.1.0"
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

            if oechem.OECalculateMolecularWeight(ligand) > 900.0:  # Units are in Dalton
                self.opt['Logger'].warn("[{}] The molecule {} seems to have a large molecular "
                                        "weight for a ligand: {:.2f} Da"
                                        .format(self.title,
                                                ligand.GetTitle(),
                                                oechem.OECalculateMolecularWeight(ligand)))

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
                self.log.info("[{}] ELF10 charge method applied to the ligand: {}".format(self.title,
                                                                                          ligand.GetTitle()))

            record.set_value(Fields.primary_molecule, ligand)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)


class LigandSetting(OERecordComputeCube):
    title = "Ligand Setting"
    version = "0.1.0"
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

    def process(self, record, port):
        try:
            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing '{}' field".format(Fields.primary_molecule.get_name()))
                raise ValueError("Missing Primary Molecule")

            ligand = record.get_value(Fields.primary_molecule)

            if oechem.OECalculateMolecularWeight(ligand) > 900.0:  # Units are in Dalton
                self.opt['Logger'].warn("[{}] The molecule {} seems to have a large molecular "
                                        "weight for a ligand: {:.2f} Da"
                                        .format(self.title,
                                                ligand.GetTitle(),
                                                oechem.OECalculateMolecularWeight(ligand)))

            record.set_value(Fields.ligand_name, ligand.GetTitle())

            ligand.SetTitle('l'+ligand.GetTitle())

            for at in ligand.GetAtoms():
                residue = oechem.OEAtomGetResidue(at)
                residue.SetName(self.args.lig_res_name)
                oechem.OEAtomSetResidue(at, residue)

            record.set_value(Fields.primary_molecule, ligand)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)
