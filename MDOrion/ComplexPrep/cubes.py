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
from datarecord import OERecord
from cuberecord import OERecordComputeCube
from cuberecord.ports import RecordInputPort

from floe.api import parameter
from openeye import oechem
from oeommtools import utils as oeommutils

from MDOrion.Standards import Fields

from floe.constants import ADVANCED


class ComplexPrepCube(OERecordComputeCube):
    title = "Complex Preparation"
    version = "0.1.0"
    classification = [["System Preparation"]]
    tags = ['Complex', 'Ligand', 'Protein']
    description = """
    This cube assembles the complex made of a protein and its docked ligands. 
    Each ligand must have just one conformer. In order to deal with multiple 
    conformers, the ligands must be processed by the “ID Setting Cube” which 
    will split ligand conformers in single conformer. In addition, each ligand 
    needs to have a ligand ID that can be set by using the “ID Setting Cube” as 
    well. The ligands must be docked to the target protein otherwise a runtime 
    error will be raised. If crystallographic water molecules are present in 
    the target protein, the water molecules that clashes with the docked ligands 
    will be removed. The ligand is identified by the ligand residue name that 
    can be set by using the cube parameter. 
    
    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the ligands
    oechem.OEDataRecord - Streamed-in of a single target protein
    

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of records with the generated complexes    
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
                                             help_text='The ligand residue name',
                                             level=ADVANCED)

    protein_port = RecordInputPort("protein_port", initializer=True)

    def begin(self):
        for record in self.protein_port:
            self.opt = vars(self.args)
            self.opt['Logger'] = self.log

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing Protein field")
                self.failure.emit(record)
                return

            protein = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                self.log.warn("Missing Protein Title field")
                self.protein_title = protein.GetTitle()[0:12]
            else:
                self.protein_title = record.get_value(Fields.title)

            self.protein = protein
            return

    def process(self, record, port):
        try:
            if port == 'intake':

                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("Missing the Primary Molecule field")

                ligand = record.get_value(Fields.primary_molecule)

                if ligand.NumConfs() > 1:
                    raise ValueError("The ligand {} has multiple conformers: {}".format(ligand.GetTitle(),
                                                                                        ligand.GetNumConfs()))

                if not record.has_value(Fields.title):
                    self.log.warn("Missing Ligand record '{}' field".format(Fields.title.get_name()))
                    ligand_title = ligand.GetTitle()[0:12]
                else:
                    ligand_title = record.get_value(Fields.title)

                if not record.has_value(Fields.id):
                    raise ValueError("Missing Ligand record '{}' field".format(Fields.id.get_name()))

                ligand_id = record.get_value(Fields.id)

                complx = self.protein.CreateCopy()
                oechem.OEAddMols(complx, ligand)

                # Split the complex in components
                protein_split, ligand_split, water, excipients = oeommutils.split(complx,
                                                                                  ligand_res_name=self.opt['lig_res_name'])

                # If the protein does not contain any atom emit a failure
                if not protein_split.NumAtoms():  # Error: protein molecule is empty
                    raise ValueError("The protein molecule does not contains atoms")

                # If the ligand does not contain any atom emit a failure
                if not ligand_split.NumAtoms():  # Error: ligand molecule is empty
                    raise ValueError("The Ligand molecule does not contains atoms")

                # Check if the ligand is inside the binding site. Cutoff distance 3A
                if not oeommutils.check_shell(ligand_split, protein_split, 3):
                    raise ValueError("The ligand is probably outside the protein binding site")

                # Removing possible clashes between the ligand and water or excipients
                if water.NumAtoms():
                    water_del = oeommutils.delete_shell(ligand, water, 1.5, in_out='in')

                if excipients.NumAtoms():
                    excipient_del = oeommutils.delete_shell(ligand, excipients, 1.5, in_out='in')

                # Reassemble the complex
                new_complex = protein_split.CreateCopy()
                oechem.OEAddMols(new_complex, ligand_split)
                if excipients.NumAtoms():
                    oechem.OEAddMols(new_complex, excipient_del)
                if water.NumAtoms():
                    oechem.OEAddMols(new_complex, water_del)

                complex_title = 'p' + self.protein_title + '_' + ligand_title

                new_complex.SetTitle(ligand.GetTitle())

                new_record = OERecord()

                # Copy all the ligand fields into the new record
                for field in record.get_fields():
                    new_record.set_value(field, record.get_value(field))

                new_record.set_value(Fields.primary_molecule, new_complex)
                new_record.set_value(Fields.title, complex_title)
                new_record.set_value(Fields.ligand, ligand)
                new_record.set_value(Fields.protein, self.protein)
                new_record.set_value(Fields.id, ligand_id)

                self.success.emit(new_record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
