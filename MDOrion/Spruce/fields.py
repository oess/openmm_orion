# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
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
from datarecord import (OEField,
                        Types)
from datarecord import OEPrimaryMolField

class SpruceFields():
    @property
    def pdb_code(self):
        return OEField('pdb_code', Types.String)

    @property
    def pdb_file(self):
        return OEField('pdb_file', Types.Blob)

    @property
    def mtz_file(self):
        return OEField('mtz_file', Types.Blob)

    @property
    def meta_data(self):
        return OEField('meta_data', Types.String)

    @property
    def pdb_url(self):
        return OEField('pdb_url', Types.String)

    @property
    def mtz_url(self):
        return OEField('mtz_url', Types.String)

    @property
    def du_vec(self):
        return OEField('du_vec', Types.BlobVec)

    @property
    def du_single(self):
        return OEField('du_single', Types.Blob)

    @property
    def du_title(self):
        return OEField('du_title', Types.String)

    @property
    def ird_vec(self):
        return OEField('ird_vec', Types.StringVec)

    @property
    def ird_single(self):
        return OEField('ird_single', Types.String)

    @property
    def receptor_vec(self):
        return OEField('receptor_vec', Types.BlobVec)

    @property
    def one_receptor(self):
        return OEField('one_receptor', Types.Chem.Mol)

    @property
    def du_ref_file(self):
        return OEField('du_ref_file', Types.Blob)

    @property
    def bu_ref_mol(self):
        return OEField('bu_ref_mol', Types.Chem.Mol)

    @property
    def ligand(self):
        return OEField('ligand', Types.Chem.Mol)

    @property
    def receptors(self):
        return OEField('receptors', Types.Chem.Mol)

    @property
    def md_complex(self):
        return OEPrimaryMolField('md_complex')