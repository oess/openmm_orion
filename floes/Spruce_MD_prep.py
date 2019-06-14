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

from floe.api import WorkFloe
from MDOrion.Spruce.prep import (PDBCodeToUrl,
                                 UrlToFile,
                                 SprucePrepAdvanced,
                                 DUtoReceptorDataset,
                                 DUtoMDDataset)

# Declare and document floe
prep_floe = WorkFloe('prep_floe', title="Spruce MD Prep Floe")

prep_floe.description = """
THIS is an ALPHA Version.

A basic floe for Spruce Protein Prep generating a protein ready for input
into MD. Specifically, residues with missing sidechain atoms have the
sidechains rebuilt, chain termini are capped with N-methyl and acetate,
and the protein is intelligently protonated.
Internal waters are retained. The cognate ligand is removed, as are
excipients, metals, and glycosylation sugars on glycosylated residues.

Known Bugs (being fixed):
- Once glycosylation sugars have been removed, residue sidechains are not yet
capped with hydrogen.
- All metals are removed even if they are important (e.g. cofactors).

Required Input Parameters:
--------------------------
PDB ID: Protein ID number

Outputs:
--------
out:  OERecord with the MD prepared protein

"""

prep_floe.classification = [["Spruce", "Examples", "Floes"]]
prep_floe.tags = ["Spruce Prep Flow", "Prep a biomolecule with Spruce in Orion."]

# Declare Cubes
pdb_url_cube = PDBCodeToUrl('pdb_url_cube')
url_download_cube = UrlToFile('url_download_cube')
prep_cube = SprucePrepAdvanced('prep_cube')
dataset_cube = DUtoReceptorDataset('dataset_cube')
dataset_md_cube = DUtoMDDataset('dataset_md_cube')

# Add cubes to floe
prep_floe.add_cubes(pdb_url_cube,
                    url_download_cube,
                    prep_cube,
                    dataset_cube,
                    dataset_md_cube)

# PDB cube parameters
pdb_url_cube.promote_parameter('pdb_code', promoted_name='pdb', title='PDB code to download to a record.', required=True)  # noqa

# pdb_url_cube.promote_parameter('pdb_codes', promoted_name='pdb_codes', title='A text file of PDB codes that will be prepped in the floe.')
# pdb_url_cube.promote_parameter('design_ref', promoted_name='design_ref', title='Reference design unit, for BU extraction, APO site extraction and alignment.')  # noqa
# pdb_url_cube.promote_parameter('biounit_ref', promoted_name='biounit_ref', title='Reference biounit for BU extraction.')  # noqa

# Spruce Prep cube parameters
# prep_cube.promote_parameter('no_protonate', promoted_name='no_protonate', title='Switch to turn off protontation')  # noqa
# prep_cube.promote_parameter('mutations', promoted_name='mutations', title='Mutations entries; format \'name,num,(insert),chain,new_name\', separate multiple by \';\'')  # noqa
# prep_cube.promote_parameter('build_sc', promoted_name='build_sc', title='Switch to turn on side chain rebuilding', default=True)  # noqa
# prep_cube.promote_parameter('cap_termini', promoted_name='cap_termini', title='Switch to turn on termini capping', default=True)  # noqa
# prep_cube.promote_parameter('no_interactions', promoted_name='no_interactions', title='Switch to not add interactions')  # noqa
# prep_cube.promote_parameter('ligand_name', promoted_name='ligand_name', title='Split only these names, separate with commas')  # noqa
# prep_cube.promote_parameter('split_cofactors', promoted_name='split_cofactors', title='Split out cofactors like ligands')  # noqa
# prep_cube.promote_parameter('cofactors', promoted_name='cofactors', title='One or more codes to be considered as cofactors')  # noqa
# prep_cube.promote_parameter('excipients', promoted_name='excipients', title='One or more codes to be considered as excipients')  # noqa
# prep_cube.promote_parameter('min_atoms', promoted_name='min_atoms', title='Min atoms for a ligand  (default: 8)')  # noqa
# prep_cube.promote_parameter('max_atoms', promoted_name='max_atoms', title='Max atoms for a ligand  (default: 100)')  # noqa
# prep_cube.promote_parameter('max_residues', promoted_name='max_residues', title='Max residues for a ligand (default: 5)')  # noqa
# prep_cube.promote_parameter('site_size', promoted_name='site_size', title='Distance between components considered in same\nsite. (small, medium, large, zero) (default: medium)')  # noqa

prep_cube.set_parameters(build_sc=True)
prep_cube.set_parameters(cap_termini=True)
prep_cube.set_parameters(mutations=None)
prep_cube.set_parameters(no_protonate=False)
prep_cube.set_parameters(no_charge=True)
prep_cube.set_parameters(no_interactions=True)
prep_cube.set_parameters(no_packing=True)
prep_cube.set_parameters(excipients=None)


# DU to Receptor cube parameters
# dataset_cube.promote_parameter('cache_grids', promoted_name='cache_grids', title='Switch to cache grids on the receptor')  # noqa
# dataset_cube.promote_parameter('keep_du', promoted_name='keep_du', title='Keep the DU on the output record.')  # noqa
# dataset_cube.promote_parameter('dataset_prefix', promoted_name='dataset_prefix', title='Prefix for the non-Orion output dataset.')  # noqa
# dataset_cube.promote_parameter('dataset_name', promoted_name='dataset_name', title='File name for the non-Orion output dataset.')  # noqa

dataset_md_cube.promote_parameter('keep_du', promoted_name='keep_du', title='Keep the DU on the output record.', default=False)  # noqa
# dataset_md_cube.promote_parameter('dataset_prefix', promoted_name='dataset_prefix', title='Prefix for the output dataset.', default='MDReady_')  # noqa
# dataset_md_cube.promote_parameter('dataset_name', promoted_name='dataset_name', title='File name for the non-Orion output dataset.', default='4WTR')  # noqa


# Connect all the ports
pdb_url_cube.success.connect(url_download_cube.intake)
url_download_cube.success.connect(prep_cube.intake)
prep_cube.success.connect(dataset_cube.intake)
prep_cube.success.connect(dataset_md_cube.intake)


if __name__ == "__main__":
    prep_floe.run()
