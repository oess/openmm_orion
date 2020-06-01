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


from openeye import oechem, oequacpac

from openmoltools.openeye import *


def assignELF10charges(molecule, max_confs=800, strictStereo=True, opt=None):
    """
     This function computes atomic partial charges for an OEMol by
     using the ELF10 method

    Parameters:
    -----------
    molecule : OEMol object
        The molecule that needs to be charged
    max_confs : integer
        The max number of conformers used to calculate the atomic partial charges
    strictStereo : bool
        a flag used to check if atoms need to have assigned stereo chemistry or not

    Return:
    -------
    mol_copy : OEMol
        a copy of the original molecule with assigned atomic partial charges
    """

    mol_copy = oechem.OEMol(molecule)

    # The passed molecule could have already conformers. If the conformer number
    # does not exceed the max_conf threshold then max_confs conformations will
    # be generated
    if not mol_copy.GetMaxConfIdx() > 200:
        # Generate up to max_confs conformers
        mol_copy = generate_conformers(mol_copy, max_confs=max_confs, strictStereo=strictStereo)

    # Assign MMFF Atom types
    if not oechem.OEMMFFAtomTypes(mol_copy):
        raise RuntimeError("MMFF atom type assignment returned errors")

    # Check for Carboxylic Acid patterns in the molecule
    smarts = '(O=)[C][O,S][H]'
    ss = oechem.OESubSearch(smarts)

    oechem.OEPrepareSearch(mol_copy, ss)
    unique_match = True

    a_match_list = []
    for match in ss.Match(mol_copy, unique_match):

        for ma in match.GetAtoms():
            a_match_list.append(ma.target)

    # Set the Carboxylic Acid torsion to zero for each generated conformers
    if a_match_list:

        if len(a_match_list) % 4 != 0:
            raise ValueError("The atom matching list must be multiple of 4")

        for i in range(0, len(a_match_list), 4):

            chunk = a_match_list[i:i + 4]

            for conf in mol_copy.GetConfs():

                conf.SetTorsion(chunk[0],
                                chunk[1],
                                chunk[2],
                                chunk[3], 0.0)

    # Try to calculate the ELF10 charges for the molecule
    quacpac_status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCELF10Charges())

    if not quacpac_status:
        opt['Logger'].warn("OEAM1BCCELF10 charge assignment failed downgrading "
                           "to OEAM1BCC charge assignment for this molecule: {}".format(mol_copy.GetTitle()))

        quacpac_status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCCharges())

    if not quacpac_status:
        raise RuntimeError("OEAssignCharges returned error code {}".format(quacpac_status))

    return mol_copy



