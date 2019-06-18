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
from floe.api import (StringParameter,
                      BooleanParameter,
                      IntegerParameter,
                      FileInputParameter)

from .fields import SpruceFields


class SpruceParameters():
    @property
    def pdb_code(self):
        return StringParameter("pdb_code",
                               title="The PDB code that will be downloaded.",
                               required=False)

    @property
    def pdb_codes(self):
        return FileInputParameter("pdb_codes",
                                  title="A text file of PDB codes (one per line) that will be downloaded.",
                                  required=False,
                                  null=True,
                                  default=None)

    @property
    def pdb_file(self):
        return FileInputParameter("pdb_file",
                                  title="The PDB file that will be used.",
                                  required=False,
                                  null=True,
                                  default=None)

    @property
    def mtz_file(self):
        return FileInputParameter("mtz_file",
                                  title="The MTZ file that will be used.",
                                  required=False,
                                  null=True,
                                  default=None)

    @property
    def design_ref(self):
        return FileInputParameter("design_ref",
                                  title='Reference design unit, for BU extraction, APO site extraction and alignment.',  # noqa
                                  required=False, default=None)

    @property
    def biounit_ref(self):
        return FileInputParameter("biounit_ref",
                                  title='Reference biounit for BU extraction.',  # noqa
                                  required=False, default=None)

    @property
    def cache_grids(self):
        return BooleanParameter("cache_grids",
                                title='Switch to cache grids on the receptor',  # noqa
                                required=True, default=False)

    @property
    def keep_du(self):
        return BooleanParameter("keep_du",
                                title="Pass the input DU along to the receptor cube",
                                required=True, default=False)

    @property
    def use_mol_field(self):
        return BooleanParameter("use_mol_field",
                                title="Output the receptor to a molecule field rather than a BlobVec.",
                                required=True, default=False)

    @property
    def no_protonate(self):
        return BooleanParameter("no_protonate",
                                title='Switch to turn off protontation',
                                required=True, default=False)

    @property
    def mutations(self):
        return StringParameter("mutations",
                                title='Mutations entries; format \'name,num,(insert),chain,new_name\', separate multiple by \';\'',  # noqa
                                required=False, default=None)

    @property
    def build_sc(self):
        return BooleanParameter("build_sc",
                                title='Switch to turn on side chain rebuilding',  # noqa
                                required=True, default=False)

    @property
    def cap_termini(self):
        return BooleanParameter("cap_termini",
                                title='Switch to turn on termini capping',
                                required=True, default=False)

    @property
    def no_interactions(self):
        return BooleanParameter("no_interactions",
                                title='Switch to not add interactions to the DU',  # noqa
                                required=True, default=False)

    @property
    def ligand_name(self):
        return StringParameter("ligand_name",
                               title='Split only these names, separate with commas',  # noqa
                               required=False, default=None)

    @property
    def split_cofactors(self):
        return BooleanParameter("split_cofactors",
                                title='Split out cofactors like ligands',  # noqa
                                required=True, default=False)

    @property
    def cofactors(self):
        return StringParameter("cofactors",
                                title='One or more codes to be considered as cofactors',  # noqa
                                required=False, default=None)

    @property
    def excipients(self):
        return StringParameter("excipients",
                               title='One or more codes to be considered as excipients',  # noqa
                               required=False, default=None)

    @property
    def min_atoms(self):
        return IntegerParameter("min_atoms",
                                title='Min atoms for a ligand  (default: 8)',
                                required=True, default=8)

    @property
    def max_atoms(self):
        return IntegerParameter("max_atoms",
                                title='Max atoms for a ligand  (default: 100)',  # noqa
                                required=True, default=100)

    @property
    def max_residues(self):
        return IntegerParameter("max_residues",
                                title='Max residues for a ligand (default: 5)',  # noqa
                                required=True, default=5)

    @property
    def max_system_atoms(self):
        return IntegerParameter("max_system_atoms",
                                title='Max number of proteins atoms allowed (default: 50000)',  # noqa
                                required=True, default=50000)

    @property
    def site_size(self):
        return StringParameter("site_size",
                               title='Distance between components considered in same\nsite. (small, medium, large, zero) (default: medium)',  # noqa
                               required=False,
                               choices=["small", "medium", "large", "zero"],
                               default='medium')
    @property
    def dataset_prefix(self):
        return StringParameter("dataset_prefix",
                               title="The prefix for the created dataset(s).",
                               required=False, default="")

    @property
    def dataset_name(self):
        return StringParameter("dataset_name",
                               title="The name of the output dataset.",
                               required=False, default="")

    @property
    def wanted_du(self):
        return StringParameter("wanted_du",
                               title="A fuzzy representation of the ligand name by which to choose the DU.",
                               required=False, default="")

    @property
    def primary_mol_field_name(self):
        return StringParameter("primary_mol_field_name",
                               title='Name of the primary molecule field ({0}, {1}) (default: {0})'.format(SpruceFields().ligand.get_name(),
                                                                                                           SpruceFields().receptors.get_name()),  # noqa
                               required=True,
                               choices=[SpruceFields().ligand.get_name(), SpruceFields().receptors.get_name()],
                               default=SpruceFields().ligand.get_name())

    @property
    def no_charge(self):
        return BooleanParameter("no_charge",
                                title='Charge the system',
                                required=False, default=True)

    @property
    def no_packing(self):
        return BooleanParameter("no_packing",
                                title='Remove or not packing residues',
                                required=False, default=True)