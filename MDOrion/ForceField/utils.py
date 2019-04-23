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


from openeye import oechem

from oeommtools import utils as oeommutils

import MDOrion

from simtk.openmm import app

import parmed

from openeye import oequacpac

from MDOrion.LigPrep import ff_utils

import numpy as np

import itertools

import os

proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS',
                   'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE',
                   'LYS', 'PHE', 'SER', 'TRP', 'VAL']

rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']


proteinff = {'Amber99SBildn': 'amber99sbildn.xml',
             'AmberFB15': 'amberfb15.xml',
             'Amber14SB': 'amber14/protein.ff14SB.xml'}

ligandff = {'Gaff': 'GAFF',
            'Gaff2': 'GAFF2',
            'Smirnoff': 'SMIRNOFF'}

solventff = {'Tip3p': 'tip3p.xml'}

otherff = {'Gaff': 'GAFF',
           'Gaff2': 'GAFF2',
           'Smirnoff': 'SMIRNOFF'}


def applyffProtein(protein, opt):
    """
    This function applies the selected force field to the
    protein

    Parameters:
    -----------
    protein: OEMol molecule
        The protein to parametrize
    opt: python dictionary
        The options used to parametrize the protein

    Return:
    -------
    protein_structure: Parmed structure instance
        The parametrized protein parmed structure
    """

    opt['Logger'].info("[{}] Protein parametrized by using: {}".format(opt['CubeTitle'],
                                                                       opt['protein_forcefield']))

    topology, positions = oeommutils.oemol_to_openmmTop(protein)

    forcefield = app.ForceField(proteinff[opt['protein_forcefield']])

    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        # Extended ff99SBildn force field
        opt['Logger'].warn("The following protein residues are not recognized by the selected FF: {} - {}"
                           "\n...Extended FF is in use".format(opt['protein_forcefield'], unmatched_residues))

        PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
        FILE_DIR = os.path.join(PACKAGE_DIR, "MDOrion/ForceField", "ffext")

        ffext_fname = os.path.join(FILE_DIR, 'amber99SBildn_ext.xml')
        forcefield = app.ForceField()
        forcefield.loadFile(ffext_fname)

        unmatched_residues = forcefield.getUnmatchedResidues(topology)

        if unmatched_residues:
            raise ValueError("Error. The following protein residues are not recognized "
                             "by the extended force field {}".format(unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False, constraints=None)
    protein_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return protein_structure


def applyffWater(water, opt):
    """
    This function applies the selected force field to the
    water

    Parameters:
    -----------
    water: OEMol molecule
        The water molecules to parametrize
    opt: python dictionary
        The options used to parametrize the water

    Return:
    -------
    water_structure: Parmed structure instance
        The parametrized water parmed structure
    """

    opt['Logger'].info("[{}] Water parametrized by using: {}".format(opt['CubeTitle'],
                                                                     opt['solvent_forcefield']))

    topology, positions = oeommutils.oemol_to_openmmTop(water)

    forcefield = app.ForceField(solventff[opt['solvent_forcefield']])

    # if opt['solvent_forcefield'] == 'tip4pew.xml':
    #     modeller = app.Modeller(topology, positions)
    #     modeller.addExtraParticles(forcefield)
    #     topology = modeller.topology
    #     positions = modeller.positions

    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        raise ValueError("The following water molecules are not recognized "
                         "by the selected force field {}: {}".format(opt['solvent_forcefield'], unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False, constraints=None)
    water_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return water_structure


def applyffExcipients(excipients, opt):
    """
    This function applies the selected force field to the
    excipients

    Parameters:
    -----------
    excipients: OEMol molecule
        The excipients molecules to parametrize
    opt: python dictionary
        The options used to parametrize the excipients

    Return:
    -------
    excipient_structure: Parmed structure instance
        The parametrized excipient parmed structure
    """

    opt['Logger'].info("[{}] Excipients parametrized by using: {} and Ions by using Amber 14".format(opt['CubeTitle'],
                                                                                                     opt['other_forcefield']))

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(excipients)


    # Ions are contained in the amberff14 force field

    exc_ff = "amber14/tip3p.xml"

    # Try to apply the selected FF on the excipients
    forcefield = app.ForceField(exc_ff)

    # List of the unrecognized excipients
    unmatched_res_list = forcefield.getUnmatchedResidues(topology)

    # Unique unrecognized excipient names
    templates = set()
    for res in unmatched_res_list:
        templates.add(res.name)

    if templates:  # Some excipients are not recognized
        opt['Logger'].warn("The following excipients are not recognized by the protein FF: {}"
                           "\nThey will be parametrized by using the FF: {}".format(templates, opt['other_forcefield']))

        # Create a bit vector mask used to split recognized from un-recognize excipients
        bv = oechem.OEBitVector(excipients.GetMaxAtomIdx())
        bv.NegateBits()

        # Dictionary containing the name and the parmed structures of the unrecognized excipients
        unrc_excipient_structures = {}

        # Dictionary used to skip already selected unrecognized excipients and count them
        unmatched_excp = {}

        # Ordered list of the unrecognized excipients
        unmatched_res_order = []

        for r_name in templates:
            unmatched_excp[r_name] = 0

        hv = oechem.OEHierView(excipients)

        for chain in hv.GetChains():
            for frag in chain.GetFragments():
                for hres in frag.GetResidues():
                    r_name = hres.GetOEResidue().GetName()
                    if r_name not in unmatched_excp:
                        continue
                    else:
                        unmatched_res_order.append(r_name)
                        if unmatched_excp[r_name]:  # Test if we have selected the unknown excipient
                            # Set Bit mask
                            atms = hres.GetAtoms()
                            for at in atms:
                                bv.SetBitOff(at.GetIdx())
                            unmatched_excp[r_name] += 1
                        else:
                            unmatched_excp[r_name] = 1
                            #  Create AtomBondSet to extract from the whole excipient system
                            #  the current selected FF unknown excipient
                            atms = hres.GetAtoms()
                            bond_set = set()
                            for at in atms:
                                bv.SetBitOff(at.GetIdx())
                                bonds = at.GetBonds()
                                for bond in bonds:
                                    bond_set.add(bond)
                            atom_bond_set = oechem.OEAtomBondSet(atms)
                            for bond in bond_set:
                                atom_bond_set.AddBond(bond)

                            # Create the unrecognized excipient OEMol
                            unrc_excp = oechem.OEMol()
                            if not oechem.OESubsetMol(unrc_excp, excipients, atom_bond_set):
                                raise ValueError("Is was not possible extract the residue: {}".format(r_name))

                            # Charge the unrecognized excipient
                            if not oequacpac.OEAssignCharges(unrc_excp, oequacpac.OEAM1BCCCharges(symmetrize=True)):
                                raise ValueError("Is was not possible to charge the extract residue: {}".format(r_name))

                            # If GAFF or GAFF2 is selected as FF check for tleap command
                            if opt['other_forcefield'] in ['GAFF', 'GAFF2']:
                                ff_utils.ParamLigStructure(oechem.OEMol(), otherff[opt['other_forcefield']]).checkTleap

                            if otherff[opt['other_forcefield']] == 'SMIRNOFF':
                                unrc_excp = oeommutils.sanitizeOEMolecule(unrc_excp)

                            # Parametrize the unrecognized excipient by using the selected FF
                            pmd = ff_utils.ParamLigStructure(unrc_excp, otherff[opt['other_forcefield']],
                                                             prefix_name=opt['prefix_name']+'_'+r_name)
                            unrc_excp_struc = pmd.parameterize()
                            unrc_excp_struc.residues[0].name = r_name
                            unrc_excipient_structures[r_name] = unrc_excp_struc

        # Recognized FF excipients
        pred_rec = oechem.OEAtomIdxSelected(bv)
        rec_excp = oechem.OEMol()
        oechem.OESubsetMol(rec_excp, excipients, pred_rec)

        if rec_excp.NumAtoms() > 0:
            top_known, pos_known = oeommutils.oemol_to_openmmTop(rec_excp)

            ff_rec = app.ForceField(exc_ff)

            try:
                omm_system = ff_rec.createSystem(top_known, rigidWater=False, constraints=None)
                rec_struc = parmed.openmm.load_topology(top_known, omm_system, xyz=pos_known)
            except:
                raise ValueError("Error in the recognised excipient parametrization")

        # Unrecognized FF excipients
        bv.NegateBits()
        pred_unrc = oechem.OEAtomIdxSelected(bv)
        unrc_excp = oechem.OEMol()
        oechem.OESubsetMol(unrc_excp, excipients, pred_unrc)

        # Unrecognized FF excipients coordinates
        oe_coord_dic = unrc_excp.GetCoords()
        unrc_coords = np.ndarray(shape=(unrc_excp.NumAtoms(), 3))
        for at_idx in oe_coord_dic:
            unrc_coords[at_idx] = oe_coord_dic[at_idx]

        # It is important the order used to assemble the structures. In order to
        # avoid mismatch between the coordinates and the structures, it is convenient
        # to use the unrecognized residue order
        unmatched_res_order_count = [(k, len(list(g))) for k, g in itertools.groupby(unmatched_res_order)]

        # Merge all the unrecognized Parmed structure
        unrc_struc = parmed.Structure()

        for pair in unmatched_res_order_count:
            res_name = pair[0]
            nums = pair[1]
            unrc_struc = unrc_struc + nums*unrc_excipient_structures[res_name]

        # Set the unrecognized coordinates
        unrc_struc.coordinates = unrc_coords

        # Set the parmed excipient structure merging
        # the unrecognized and recognized parmed
        # structures together
        if rec_excp.NumAtoms() > 0:
            excipients_structure = unrc_struc + rec_struc
        else:
            excipients_structure = unrc_struc

        # return excipients_structure
    else:  # All the excipients are recognized by the selected FF
        omm_system = forcefield.createSystem(topology, rigidWater=False, constraints=None)
        excipients_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

        # return excipients_structure

    return excipients_structure


def applyffLigand(ligand, opt):
    """
    This function applies the selected force field to the
    ligand

    Parameters:
    -----------
    ligand: OEMol molecule
        The ligand molecule to parametrize
    opt: python dictionary
        The options used to parametrize the ligand

    Return:
    -------
    ligand_structure: Parmed structure instance
        The parametrized ligand parmed structure
    """

    # Check TLeap
    if opt['ligand_forcefield'] in ['GAFF', 'GAFF2']:
        ff_utils.ParamLigStructure(oechem.OEMol(), ligandff[opt['ligand_forcefield']]).checkTleap

    # Parametrize the Ligand
    pmd = ff_utils.ParamLigStructure(ligand, ligandff[opt['ligand_forcefield']], prefix_name=opt['prefix_name'])
    ligand_structure = pmd.parameterize()
    ligand_structure.residues[0].name = opt['lig_res_name']
    opt['Logger'].info("[{}] Ligand parametrized by using: {}".format(opt['CubeTitle'],
                                                                      opt['ligand_forcefield']))

    return ligand_structure
