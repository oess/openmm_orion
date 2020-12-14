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


import traceback
from datarecord import OERecord

from oeommtools import utils as oeommutils

from MDOrion.Standards import Fields

from floe.api import ComputeCube

from orionplatform.mixins import RecordPortsMixin
from orionplatform.ports import RecordInputPort

from openeye import oechem


class ComplexPrepCube(RecordPortsMixin, ComputeCube):
    title = "Complex Preparation"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['Complex', 'Ligand', 'Protein']
    description = """
    This Cube assembles the complex made of a protein and its docked ligands. 
    Each ligand must have just one conformer. In order to deal with multiple 
    conformers, the ligands must be processed by the “ID Setting Cube” which 
    will split ligand conformers in single conformer. In addition, each ligand 
    needs to have a ligand ID that can be set by using the “ID Setting Cube” as 
    well. The ligands must be docked to the target protein otherwise a runtime 
    error will be raised. If crystallographic water molecules are present in 
    the target protein, the water molecules that clashes with the docked ligands 
    will be removed. The ligand is identified by the ligand residue name that 
    can be set by using the cube parameter. 
    """

    uuid = "be2ac138-22ae-4412-9c38-886472c496b9"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_port = RecordInputPort("protein_port", initializer=True)

    def begin(self):
        for record in self.protein_port:

            self.opt = vars(self.args)
            self.opt['Logger'] = self.log

            if not record.has_value(Fields.md_components):
                raise ValueError("MD Components Field is missing")

            self.md_components = record.get_value(Fields.md_components)

        return

    def process(self, record, port):
        try:
            if port == 'intake':

                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("Missing the ligand primary molecule field")

                ligand = record.get_value(Fields.primary_molecule)

                if ligand.NumConfs() > 1:
                    raise ValueError("The ligand {} has multiple conformers: {}".format(ligand.GetTitle(),
                                                                                        ligand.GetNumConfs()))

                if not record.has_value(Fields.title):
                    self.log.warn("Missing title field '{}' field; improvising".format(Fields.title.get_name()))
                    ligand_title = ligand.GetTitle()[0:12]
                else:
                    ligand_title = record.get_value(Fields.title)

                protein = self.md_components.get_protein

                self.md_components.set_ligand(ligand)

                # Check if the ligand is inside the binding site. Cutoff distance 3A
                if not oeommutils.check_shell(ligand, protein, 3):
                    raise ValueError("The Ligand is probably outside the Protein binding site")

                # Remove Steric Clashes between the ligand and the other System components
                for comp_name, comp in self.md_components.get_components.items():

                    # Skip clashes between the ligand itself and the protein
                    if comp_name in ['ligand', 'protein']:
                        continue

                    # Remove Metal clashes if the distance between the metal and the ligand
                    # is less than 1A
                    elif comp_name == 'metals':
                        metal_del = oeommutils.delete_shell(ligand, comp, 1.0, in_out='in')

                        if metal_del.NumAtoms() != comp.NumAtoms():
                            self.opt['Logger'].info(
                                "Detected steric-clashes between the ligand {} and metals".format(ligand_title))

                            self.md_components.set_metals(metal_del)
                            # Remove  clashes if the distance between the selected component and the ligand
                            # is less than 1.5A
                    else:
                        comp_del = oeommutils.delete_shell(ligand, comp, 1.5, in_out='in')

                        if comp_del.NumAtoms() != comp.NumAtoms():
                            self.opt['Logger'].info(
                                "Detected steric-clashes between the ligand {} and component {}".format(
                                    ligand_title,
                                    comp_name))

                            self.md_components.set_component_by_name(comp_name, comp_del)

                complex_title = 'p' + self.md_components.get_title + '_l' + ligand_title

                mdcomp = self.md_components.copy
                mdcomp.set_title(complex_title)

                # Check Ligand
                lig_check = mdcomp.get_ligand
                smi_lig_check = oechem.OECreateSmiString(lig_check)
                smi_ligand = oechem.OECreateSmiString(ligand)

                if smi_ligand != smi_lig_check:
                    raise ValueError("Ligand IsoSmiles String check failure: {} vs {}".format(smi_lig_check, smi_ligand))

                # the ligand is the primary molecule
                new_record = OERecord(record)

                new_record.set_value(Fields.title, complex_title)
                new_record.set_value(Fields.ligand, ligand)
                new_record.set_value(Fields.protein, protein)

                # Check Protein Name
                if protein.GetTitle():
                    protein_name = protein.GetTitle()
                else:
                    protein_name = "prot"

                new_record.set_value(Fields.protein_name, protein_name)
                new_record.set_value(Fields.md_components, mdcomp)

                self.success.emit(new_record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
