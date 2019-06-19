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
from MDOrion.Spruce.prep import (SprucePrepBasic,
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
PDB File: Protein PDB file
MTZ File: Protein MTZ file

Outputs:
--------
out:  OERecord with the MD prepared protein
      OERecord with the OE Receptor
      OERecord with the Protein Bio-Unit

"""

prep_floe.classification = [["Spruce", "Examples", "Floes"]]
prep_floe.tags = ["Spruce Prep Flow", "Prep a biomolecule with Spruce in Orion."]

# Declare Cubes
prep_cube = SprucePrepBasic('prep_cube')

dataset_cube = DUtoReceptorDataset('dataset_cube')
dataset_md_cube = DUtoMDDataset('dataset_md_cube')

# Add cubes to floe
prep_floe.add_cubes(prep_cube,
                    dataset_cube,
                    dataset_md_cube)

# PDB cube parameters
prep_cube.promote_parameter('pdb_code', promoted_name='pdb_code', title='The protein PDB code')
prep_cube.promote_parameter('pdb_file', promoted_name='pdb_file', title='The protein PDB file')
prep_cube.promote_parameter('mtz_file', promoted_name='mtz_file', title='The protein MTZ file')


prep_cube.set_parameters(build_sc=True)
prep_cube.set_parameters(cap_termini=True)
prep_cube.set_parameters(mutations=None)
prep_cube.set_parameters(no_protonate=False)
prep_cube.set_parameters(no_charge=True)
prep_cube.set_parameters(no_interactions=True)
prep_cube.set_parameters(no_packing=True)
prep_cube.set_parameters(excipients=None)

dataset_md_cube.promote_parameter('keep_du', promoted_name='keep_du', title='Keep the DU on the output record.', default=False)

# Connect all the ports
prep_cube.success.connect(dataset_cube.intake)
prep_cube.success.connect(dataset_md_cube.intake)


if __name__ == "__main__":
    prep_floe.run()
