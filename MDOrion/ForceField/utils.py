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
             'Amber99SB': 'amber99sb.xml',
             'Amber14SB': 'amber14/protein.ff14SB.xml',
             'AmberFB15': 'amberfb15.xml'}

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
    numparts, partlist = oechem.OEDetermineComponents(excipients)
    pred = oechem.OEPartPredAtom(partlist)

    part_mols = []
    for i in range(1, numparts + 1):
        pred.SelectPart(i)
        partmol = oechem.OEMol()
        oechem.OESubsetMol(partmol, excipients, pred)
        part_mols.append(partmol)

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(excipients)

    map_omm_to_oe = {omm_res: oe_mol for omm_res, oe_mol in zip(topology.residues(), part_mols)}

    # Ions are contained in the amberff14 force field
    exc_ff = "amber14/tip3p.xml"

    # Try to apply the selected FF on the excipients
    forcefield = app.ForceField(exc_ff)
    #
    # List of the unrecognized excipients
    unmatched_res_list = forcefield.getUnmatchedResidues(topology)

    matched_res_list = []
    for res in topology.residues():
        if res in unmatched_res_list:
            continue
        else:
            matched_res_list.append(res)

    # Map OpenMM residue to their parmed structure
    map_omm_res_to_pmd = dict()

    # Unique residue templates
    map_template_to_pmd = dict()

    # Matched Residue Parametrization
    bondedToAtom = forcefield._buildBondedToAtomList(topology)
    for res in matched_res_list:
        template, matches = forcefield._getResidueTemplateMatches(res, bondedToAtom)
        if template in map_template_to_pmd:
            map_omm_res_to_pmd[res] = map_template_to_pmd[template]
        else:
            res_top, res_pos = oeommutils.oemol_to_openmmTop(map_omm_to_oe[res])
            try:
                res_omm_system = forcefield.createSystem(res_top, rigidWater=False, constraints=None)
                res_pmd = parmed.openmm.load_topology(res_top, res_omm_system, xyz=res_pos)
            except:
                raise ValueError("Error in the recognised excipient residue parametrization {}".format(res))

            map_template_to_pmd[template] = res_pmd
            map_omm_res_to_pmd[res] = res_pmd

    # print(map_omm_res_to_pmd)

    # UnMatched Residue Parametrization
    for ures in unmatched_res_list:

        oe_mol = map_omm_to_oe[ures]

        # Charge the unrecognized excipient
        if not oequacpac.OEAssignCharges(oe_mol, oequacpac.OEAM1BCCCharges(symmetrize=True)):
            raise ValueError("Is was not possible to charge the extract residue: {}".format(ures))

        # If GAFF or GAFF2 is selected as FF check for tleap command
        # This is check for the tleap command only
        if opt['other_forcefield'] in ['GAFF', 'GAFF2']:
            ff_utils.ParamLigStructure(oechem.OEMol(), otherff[opt['other_forcefield']]).checkTleap

        if otherff[opt['other_forcefield']] == 'SMIRNOFF':
            oe_mol = oeommutils.sanitizeOEMolecule(oe_mol)

        # Parametrize the unrecognized excipient by using the selected FF
        pmd = ff_utils.ParamLigStructure(oe_mol, otherff[opt['other_forcefield']],
                                         prefix_name=opt['prefix_name'] + '_' + ures.name)

        ures_pmd = pmd.parameterize()
        map_omm_res_to_pmd[ures] = ures_pmd

    # Excipient Parmed Structure
    excipients_pmd = parmed.Structure()
    for res in topology.residues():
        excipients_pmd += map_omm_res_to_pmd[res]

    if len(excipients_pmd.atoms) != excipients.NumAtoms():
        raise ValueError(
            "Excipient OE molecule and Excipient Parmed structure number of atoms mismatch {} vs {}".format(
                len(excipients_pmd.atoms), excipients.NumAtoms()))

    # Set the positions
    oe_exc_coord_dic = excipients.GetCoords()
    exc_coords = np.ndarray(shape=(excipients.NumAtoms(), 3))
    for at_idx in oe_exc_coord_dic:
        exc_coords[at_idx] = oe_exc_coord_dic[at_idx]

    excipients_pmd.coordinates = exc_coords

    return excipients_pmd


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


def clean_tags(molecule):
    """
    This function remove tags that could cause problems along the MD Analysis stage.
    In particular Hint interactions, Style and PDB data are removed.

    Parameters:
    -----------
    molecule: OEMol molecule
        The molecule to clean

    Return:
    -------
    molecule: OEMol molecule
        The cleaned molecule
    """

    oechem.OEDeleteInteractionsHintSerializationData(molecule)
    oechem.OEDeleteInteractionsHintSerializationIds(molecule)
    oechem.OEClearStyle(molecule)
    oechem.OEClearPDBData(molecule)

    return molecule
