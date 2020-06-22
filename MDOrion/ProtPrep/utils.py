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


from openeye import oechem


def ss_bond_fix(protein_arg):

    protein = oechem.OEMol(protein_arg)

    nn = oechem.OENearestNbrs(protein, 2.5)

    disulfide_atom_list = []

    for atom in protein.GetAtoms(oechem.OEIsSulfur()):
        for nbrs in nn.GetNbrs(atom):
            if atom.GetIdx() >= nbrs.GetBgn().GetIdx():
                continue
            if nbrs.GetBgn().GetAtomicNum() != oechem.OEElemNo_S:
                continue

            dis_set = {atom.GetIdx(), nbrs.GetBgn().GetIdx()}

            disulfide_atom_list.append(dis_set)

            # print(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()), atom.GetIdx(),
            #       oechem.OEGetAtomicSymbol(nbrs.GetBgn().GetAtomicNum()), nbrs.GetBgn().GetIdx())

    bond_list_set = [{bond.GetBgn().GetIdx(), bond.GetEnd().GetIdx()} for bond in protein.GetBonds()]

    missing_dis = []

    for s in disulfide_atom_list:
        if s not in bond_list_set:
            missing_dis.append(s)

    # print(missing_dis)

    for s in missing_dis:
        tmp_list = list(s)
        oe_at0 = protein.GetAtom(oechem.OEHasAtomIdx(tmp_list[0]))
        oe_at1 = protein.GetAtom(oechem.OEHasAtomIdx(tmp_list[1]))

        protein.NewBond(oe_at0, oe_at1, 1)

        delete_h_list = []
        for nbor in oe_at0.GetAtoms(oechem.OEIsHydrogen()):
            delete_h_list.append(nbor)

        for nbor in oe_at1.GetAtoms(oechem.OEIsHydrogen()):
            delete_h_list.append(nbor)

        for atd in delete_h_list:
            protein.DeleteAtom(atd)

        oe_at0.SetImplicitHCount(0)
        oe_at1.SetImplicitHCount(0)

    return protein

