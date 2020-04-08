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

from openeye import oechem
from oeommtools import utils as oeommutils

from MDOrion.Standards import Fields

from floe.constants import ADVANCED

from floe.api import (parameters,
                      ComputeCube)

from orionplatform.mixins import RecordPortsMixin
from orionplatform.ports import RecordInputPort

from openeye import oespruce


class ComplexPrepCube(RecordPortsMixin, ComputeCube):
    title = "Complex Preparation"
    # version = "0.1.4"
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

            if record.has_value(Fields.design_unit):
                du = record.get_value(Fields.design_unit)
                protein_flask = oechem.OEMol()
                if not self.du.GetProtein(protein_flask):
                    raise ValueError("The protein cannot be extracted from the Design Unit")
                self.du = du
            else:

                if not record.has_value(Fields.flask):
                    self.log.error("Missing Protein flask field")
                    self.failure.emit(record)
                    return

                protein_flask = record.get_value(Fields.flask)

            if not record.has_value(Fields.title):
                self.log.warn("Missing Protein Title field")
                self.protein_flask_title = protein_flask.GetTitle()[0:12]
            else:
                self.protein_flask_title = record.get_value(Fields.title)

            self.protein_flask = protein_flask
            return

    def process(self, record, port):
        try:
            if port == 'intake':

                if not record.has_value(Fields.flask):
                    raise ValueError("Missing the ligand flask field")

                ligand = record.get_value(Fields.flask)

                if ligand.NumConfs() > 1:
                    raise ValueError("The ligand {} has multiple conformers: {}".format(ligand.GetTitle(),
                                                                                        ligand.GetNumConfs()))

                if not record.has_value(Fields.title):
                    self.log.warn("Missing title field '{}' field; improvising".format(Fields.title.get_name()))
                    ligand_title = ligand.GetTitle()[0:12]
                else:
                    ligand_title = record.get_value(Fields.title)

                # Create a Design Unit
                if not record.has_value(Fields.design_unit):

                    du = oechem.OEDesignUnit()

                    build_opts = oespruce.OEDesignUnitBuildOptions()
                    build_opts.SetBuildSidechains(False)
                    build_opts.SetBuildLoops(False)
                    build_opts.SetCapCTermini(False)
                    build_opts.SetCapNTermini(False)
                    build_opts.SetDeleteClashingSolvent(False)

                    enum_opts = oespruce.OEDesignUnitEnumerateSitesOptions()
                    enum_opts.SetAddInteractionHints(False)
                    enum_opts.SetAddStyle(False)

                    prep_opts = oespruce.OEDesignUnitPrepOptions()
                    prep_opts.SetProtonate(False)
                    prep_opts.SetAssignPartialChargesAndRadii(False)
                    prep_opts.SetBuildOptions(build_opts)
                    prep_opts.SetEnumerateSitesOptions(enum_opts)

                    split_opts = oespruce.OEDesignUnitSplitOptions()
                    split_opts.SetMakePackingResidues(False)

                    bio_opts = oespruce.OEBioUnitExtractionOptions()
                    bio_opts.SetSuperpose(False)

                    du_opts = oespruce.OEMakeDesignUnitOptions(split_opts, prep_opts, bio_opts)

                    if not oespruce.OEMakeDesignUnit(du, self.protein_flask, ligand, du_opts):
                        raise ValueError("It was not possible to generate the Design Unit")

                    flask_from_du = oechem.OEMol()

                    if not du.GetComponents(flask_from_du, oechem.OEDesignUnitComponents_All ^ oechem.OEDesignUnitComponents_Ligand):
                        raise ValueError("Recovering all the Design Unit Components Failed")

                    if flask_from_du.NumAtoms() != self.protein_flask.NumAtoms():
                        raise ValueError("Flask Design Unit and User protein flask number "
                                         "of atoms mismatch: {} vs {}".format(flask_from_du.NumAtoms(),
                                                                              self.protein_flask.NumAtoms()))
                    self.du = du

                else:  # Update the Design Unit Ligand Component
                    if not oechem.OEUpdateDesignUnit(self.du, ligand, oechem.OEDesignUnitComponents_Ligand):
                        raise ValueError("Could not add the ligand to the Design Unit")

                protein_comp = oechem.OEMol()
                ligand_comp = oechem.OEMol()

                # Ligand and Protein component Check
                if not self.du.GetProtein(protein_comp):
                    raise ValueError("The protein cannot be extracted from the Design Unit")

                if not self.du.GetLigand(ligand_comp):
                    raise ValueError("The ligand cannot be extracted from the Design Unit")

                # If the protein does not contain any atom emit a failure
                if not protein_comp.NumAtoms():  # Error: protein molecule is empty
                    raise ValueError("The Protein component does not contains atoms")

                # If the ligand does not contain any atom emit a failure
                if not ligand_comp.NumAtoms():  # Error: ligand molecule is empty
                    raise ValueError("The Ligand component does not contains atoms")

                # Check if the ligand is inside the binding site. Cutoff distance 3A
                if not oeommutils.check_shell(ligand_comp, protein_comp, 3):
                    raise ValueError("The Ligand is probably outside the Protein binding site")

                # Remove Steric Clashes between the ligand and the other System components
                for pair in self.du.GetTaggedComponents():

                    comp_name = pair[0]
                    comp_id = self.du.GetComponentID(comp_name)
                    comp = pair[1]

                    # Skip clashes between the ligand itself and the protein
                    if comp_id == oechem.OEDesignUnitComponents_Protein \
                            or comp_id == oechem.OEDesignUnitComponents_Ligand:
                        continue

                    # Remove Metal clashes if the distance between the metal and the ligand
                    # is less than 1A
                    elif comp_id == oechem.OEDesignUnitComponents_Metals:
                        metal_del = oeommutils.delete_shell(ligand_comp, comp, 1.0, in_out='in')

                        if metal_del.NumAtoms() != comp.NumAtoms():
                            self.opt['Logger'].info(
                                "Detected steric-clashes between the ligand {} and metals".format(ligand_title))

                        if not oechem.OEUpdateDesignUnit(self.du, metal_del,
                                                         oechem.OEDesignUnitComponents_Metals):
                            raise ValueError("Could not add the un-clashed Metals to the Design Unit")

                    # Remove  clashes if the distance between the selected component and the ligand
                    # is less than 1.5A
                    else:
                        comp_del = oeommutils.delete_shell(ligand_comp, comp, 1.5, in_out='in')

                        if comp_del.NumAtoms() != comp.NumAtoms():
                            self.opt['Logger'].info(
                                "Detected steric-clashes between the ligand {} and component {}".format(
                                    ligand_title,
                                    comp_name))
                        if not oechem.OEUpdateDesignUnit(self.du, comp_del, comp_id):
                            raise ValueError("Could not add the un-clashed component to the Design Unit")

                # TODO check parametrization at this stage on the DU components

                flask = oechem.OEMol()

                if not self.du.GetComponents(flask, oechem.OEDesignUnitComponents_All):
                    raise ValueError("Recovering all the Design Unit Components Failed")

                complex_title = 'p' + self.protein_flask_title + '_l' + ligand_title
                self.du.SetTitle(complex_title)

                flask.SetTitle(complex_title)

                # the ligand is the primary molecule
                new_record = OERecord(record)

                new_record.set_value(Fields.flask, flask)
                new_record.set_value(Fields.title, complex_title)
                new_record.set_value(Fields.ligand, ligand)
                new_record.set_value(Fields.protein, protein_comp)
                new_record.set_value(Fields.protein_name, self.protein_flask_title)
                new_record.set_value(Fields.design_unit, self.du)

                self.success.emit(new_record)

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
