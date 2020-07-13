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

from oeommtools import utils as oeommutils

import numpy as np

import tempfile

import parmed

import openmoltools

# from openmoltools.openeye import *

import os, shutil

from simtk.openmm.app import AmberInpcrdFile, AmberPrmtopFile
from simtk.openmm import app

from openforcefield.typing.engines.smirnoff import ForceField

from pkg_resources import resource_filename

from openforcefield.topology import Topology, Molecule

from MDOrion.LigPrep.ff_utils import assignELF10charges

from MDOrion.ForceField.ff_library import ff_library

import mdtraj.utils

import subprocess

# TODO TEMPORARY SOLUTION FOR OPENMOLTOOLS BUG
# https://github.com/choderalab/openmoltools/issues/299


def run_tleap(molecule_name,
              gaff_mol2_filename,
              frcmod_filename,
              prmtop_filename=None,
              inpcrd_filename=None,
              leaprc='leaprc.gaff'):

    """Run AmberTools tleap to create simulation files for AMBER

    Parameters
    ----------
    molecule_name : str
        The name of the molecule
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    prmtop_filename : str, optional, default=None
        Amber prmtop file produced by tleap, defaults to molecule_name
    inpcrd_filename : str, optional, default=None
        Amber inpcrd file produced by tleap, defaults to molecule_name
    leaprc : str, optional, default = 'leaprc.gaff'
        Optionally, specify alternate leaprc to use, such as `leaprc.gaff2`

    Returns
    -------
    prmtop_filename : str
        Amber prmtop file produced by tleap
    inpcrd_filename : str
        Amber inpcrd file produced by tleap
    """
    if prmtop_filename is None:
        prmtop_filename = "%s.prmtop" % molecule_name
    if inpcrd_filename is None:
        inpcrd_filename = "%s.inpcrd" % molecule_name

    # Get absolute paths for input/output
    gaff_mol2_filename = os.path.abspath(gaff_mol2_filename)
    frcmod_filename = os.path.abspath(frcmod_filename)
    prmtop_filename = os.path.abspath(prmtop_filename)
    inpcrd_filename = os.path.abspath(inpcrd_filename)

    # Work in a temporary directory, on hard coded filenames,
    # to avoid any issues AMBER may have with spaces and other special characters in filenames
    with mdtraj.utils.enter_temp_directory():
        shutil.copy(gaff_mol2_filename, 'file.mol2')
        shutil.copy(frcmod_filename, 'file.frcmod')

        tleap_input = """
    source oldff/leaprc.ff99SB
    source %s
    LIG = loadmol2 file.mol2
    loadamberparams file.frcmod
    check LIG
    saveamberparm LIG out.prmtop out.inpcrd
    quit

""" % leaprc

        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_input)
        file_handle.close()

        cmd = "tleap -f %s " % file_handle.name

        subprocess.getoutput(cmd)

        # Copy back target files
        shutil.copy('out.prmtop', prmtop_filename)
        shutil.copy('out.inpcrd', inpcrd_filename)

    return prmtop_filename, inpcrd_filename


class ParamMolStructure(object):
    """
    Generates parametrized ParmEd structure of the molecule with a chosen force field

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The openeye molecule to be parameterized
    forcefield : str
        String specifying the forcefield parameters to be used
    prefix_name : str
        String specifying the output prefix filename

    Returns
    ---------
    packedmol : openeye.oechem.OEMol
        Openeye molecule with the ParmEd Structure attached.
    """

    def __init__(self, molecule, forcefield, prefix_name='ligand', delete_out_files=True, recharge=False):
        if forcefield not in list(ff_library.ligandff.values()):
            raise RuntimeError('The selected ligand force field is not '
                               'supported {}. Available {}'.format(forcefield, list(ff_library.ligandff.keys())))
        else:
            self.molecule = molecule
            self.forcefield = str(forcefield).strip()
            self.structure = None
            self.prefix_name = prefix_name
            self.delete_out_files = delete_out_files
            self.recharge = recharge

    def checkCharges(self, molecule):
        # Check that molecule is charged.
        is_charged = False
        for atom in molecule.GetAtoms():
            if atom.GetPartialCharge() != 0.0:
                is_charged = True

        if is_charged and self.recharge:
            is_charged = False

        return is_charged

    def getSmirnoffStructure(self, molecule=None):
        if not molecule:
            molecule = self.molecule

        if not self.checkCharges(molecule):
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")
            molecule = assignELF10charges(molecule)

        if self.forcefield == ff_library.ligandff['Smirnoff99Frosst']:

            fffn = resource_filename('openforcefield', os.path.join('data', 'test_forcefields/' + self.forcefield))

            if not os.path.exists(fffn):
                raise ValueError(
                    "Sorry! {} does not exist. If you just added it, you'll have to re-install".format(fffn))

            with open(fffn) as ffxml:
                ff = ForceField(ffxml, allow_cosmetic_attributes=True)

        elif self.forcefield in [ff_library.ligandff['OpenFF_1.0.0'],
                                 ff_library.ligandff['OpenFF_1.1.1'],
                                 ff_library.ligandff['OpenFF_1.2.0']]:

            ff = ForceField(self.forcefield, allow_cosmetic_attributes=True)

        else:
            raise ValueError("Force Field not Supported: {}".format(self.forcefield))

        mol_off = Molecule.from_openeye(molecule, allow_undefined_stereo=True)
        topology = Topology.from_molecules([mol_off])

        omm_sys = ff.create_openmm_system(topology, charge_from_molecules=[mol_off])

        # omm_top = generateTopologyFromOEMol(molecule)
        # positions = mol_off.conformers[0]

        omm_top, positions = oeommutils.oemol_to_openmmTop(molecule)

        pmd_structure = parmed.openmm.load_topology(omm_top, omm_sys, xyz=positions)

        return pmd_structure

    def getGaffStructure(self, molecule=None, forcefield=None):
        if not molecule:
            molecule = self.molecule

        if not forcefield:
            forcefield = self.forcefield

        if not molecule:
            molecule = self.molecule

        if not self.checkCharges(molecule):
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")
            molecule = assignELF10charges(molecule)

        # Write out mol to a mol2 file to process via AmberTools
        mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
        mol2filename = mol2file.name

        with oechem.oemolostream(mol2filename) as ofs:
            oechem.OEWriteConstMolecule(ofs, molecule)

        prefix = self.prefix_name + '_' + os.path.basename(mol2filename).split('.')[0]

        gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber(prefix, mol2filename,
                                                                                 gaff_version=forcefield.lower(),
                                                                                 charge_method=None)
        # Run tleap using specified forcefield
        # prmtop, inpcrd = openmoltools.amber.run_tleap(self.prefix_name, gaff_mol2_filename,
        #                                               frcmod_filename,
        #                                               leaprc='leaprc.{}'.format(forcefield.lower()))
        #
        prmtop, inpcrd = run_tleap(prefix,
                                   gaff_mol2_filename,
                                   frcmod_filename,
                                   leaprc='leaprc.{}'.format(forcefield.lower()))

        # TODO Load via ParmEd: This is causing Problems
        #  Merging two structures (OpenMM PMD structure and
        #  Amber PMD Structure): The NB exception list is messed up
        # molecule_structure = parmed.amber.AmberParm(prmtop, inpcrd)

        # TODO MODIFIED BY GAC
        omm_prmtop = AmberPrmtopFile(prmtop)
        omm_inpcrd = AmberInpcrdFile(inpcrd)

        omm_system = omm_prmtop.createSystem(nonbondedMethod=app.NoCutoff)

        molecule_structure = parmed.openmm.load_topology(omm_prmtop.topology, omm_system, xyz=omm_inpcrd.positions)

        if self.delete_out_files:
            os.remove(gaff_mol2_filename)
            os.remove(frcmod_filename)
            os.remove(prmtop)
            os.remove(inpcrd)

        return molecule_structure

    def parameterize(self):

        if self.forcefield in [ff_library.ligandff['OpenFF_1.0.0'],
                               ff_library.ligandff['OpenFF_1.1.1'],
                               ff_library.ligandff['OpenFF_1.2.0'],
                               ff_library.ligandff['Smirnoff99Frosst']]:

            structure = self.getSmirnoffStructure()

        elif self.forcefield in ['GAFF', 'GAFF2']:
            structure = self.getGaffStructure()

        self.structure = structure
        return self.structure


def parametrize_component(component, component_ff, other_ff):

    component_copy = oechem.OEMol(component)

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(component_copy)

    # Try to apply the selected FF on the component
    if isinstance(component_ff, list):
        forcefield = app.ForceField()
        for f in component_ff:
            forcefield.loadFile(f)
    else:
        forcefield = app.ForceField(component_ff)

    # List of the unrecognized component
    unmatched_res_list = forcefield.getUnmatchedResidues(topology)

    # Try to parametrize the whole flask
    if not unmatched_res_list:
        omm_components = forcefield.createSystem(topology, rigidWater=False, constraints=None)
        components_pmd = parmed.openmm.load_topology(topology, omm_components, xyz=positions)
        return components_pmd

    # Extract The different non bonded parts
    numparts, parts = oechem.OEDetermineComponents(component_copy)
    pred = oechem.OEPartPredAtom(parts)

    part_mols = dict()
    for i in range(1, numparts + 1):
        pred.SelectPart(i)
        partmol = oechem.OEMol()
        oechem.OESubsetMol(partmol, component_copy, pred)
        part_mols[i] = partmol

    part = -1
    part_numbers = []
    for at in component_copy.GetAtoms():
        if parts[at.GetIdx()] != part:
            part = parts[at.GetIdx()]
            part_numbers.append(part)
        # print("atom %d is in part %d" % (at.GetIdx(), parts[at.GetIdx()]))

    part_numbers_resolved = dict()
    for i in range(0, len(part_numbers)):
        if part_numbers[i] in part_numbers_resolved:
            continue
        else:
            part_numbers_resolved[part_numbers[i]] = part_numbers[i]
            if i == len(part_numbers) - 1:
                break
            for j in range(i+1, len(part_numbers)):
                if part_numbers[j] in part_numbers_resolved:
                    continue
                else:
                    smiles_i = oechem.OECreateCanSmiString(part_mols[part_numbers[i]])
                    smiles_j = oechem.OECreateCanSmiString(part_mols[part_numbers[j]])

                    if smiles_i == smiles_j:
                        part_numbers_resolved[part_numbers[j]] = part_numbers[i]
                    else:
                        continue

    ordered_parts = [part_numbers_resolved[pn] for pn in part_numbers]

    # print(ordered_parts)

    unique_parts = set(ordered_parts)

    part_to_pmd = dict()
    for part_i in unique_parts:

        mol_part_i = part_mols[part_i]
        omm_top_part_i, omm_pos_part_i = oeommutils.oemol_to_openmmTop(mol_part_i)

        if not forcefield.getUnmatchedResidues(omm_top_part_i):

            omm_sys_part_i = forcefield.createSystem(omm_top_part_i, rigidWater=False, constraints=None)

            pmd_part_i = parmed.openmm.load_topology(omm_top_part_i, omm_sys_part_i, xyz=omm_pos_part_i)
            # print("Parametrize full part: {}".format(part_i))
        else:
            mol_part_i_clean = oeommutils.sanitizeOEMolecule(mol_part_i)
            pmd_part_i = ParamMolStructure(mol_part_i_clean, other_ff,
                                           prefix_name="Part" + '_' + str(part_i),
                                           recharge=True).parameterize()
            # print("Parametrize part: {}".format(part_i))
        part_to_pmd[part_i] = pmd_part_i

    # Component Parmed Structure
    component_pmd = parmed.Structure()
    for part_i in ordered_parts:
        component_pmd += part_to_pmd[part_i]

    # Checking
    for parm_at, oe_at in zip(component_pmd.atoms, component_copy.GetAtoms()):

        if parm_at.atomic_number != oe_at.GetAtomicNum():
            raise ValueError(
                "Atomic number mismatch between the Parmed and the OpenEye topologies: {} - {}".
                format(parm_at.atomic_number, oe_at.GetAtomicNum()))

    # Set the positions
    oe_comp_coord_dic = component_copy.GetCoords()
    comp_coords = np.ndarray(shape=(component_copy.NumAtoms(), 3))
    for at_idx in oe_comp_coord_dic:
        comp_coords[at_idx] = oe_comp_coord_dic[at_idx]

    component_pmd.coordinates = comp_coords

    return component_pmd


# def parametrize_component(component, component_ff, other_ff):
#
#     component_copy = oechem.OEMol(component)
#
#     # OpenMM topology and positions from OEMol
#     topology, positions = oeommutils.oemol_to_openmmTop(component_copy)
#
#     # Try to apply the selected FF on the component
#     if isinstance(component_ff, list):
#         forcefield = app.ForceField()
#         for f in component_ff:
#             forcefield.loadFile(f)
#     else:
#         forcefield = app.ForceField(component_ff)
#
#     # List of the unrecognized component
#     unmatched_res_list = forcefield.getUnmatchedResidues(topology)
#
#     if not unmatched_res_list:
#         omm_components = forcefield.createSystem(topology, rigidWater=False, constraints=None)
#         components_pmd = parmed.openmm.load_topology(topology, omm_components, xyz=positions)
#         return components_pmd
#
#     numparts, partlist = oechem.OEDetermineComponents(component_copy)
#     pred = oechem.OEPartPredAtom(partlist)
#
#     part_mols = []
#     for i in range(1, numparts + 1):
#         pred.SelectPart(i)
#         partmol = oechem.OEMol()
#         oechem.OESubsetMol(partmol, component_copy, pred)
#         part_mols.append(partmol)
#
#     map_omm_to_oe = {omm_res: oe_mol for omm_res, oe_mol in zip(topology.residues(), part_mols)}
#
#     matched_res_list = []
#     for res in topology.residues():
#         if res in unmatched_res_list:
#             continue
#         else:
#             matched_res_list.append(res)
#
#     # Map OpenMM residue to their parmed structure
#     map_omm_res_to_pmd = dict()
#
#     # Unique residue templates
#     map_template_to_pmd = dict()
#
#     # Matched Residue Parametrization
#     bondedToAtom = forcefield._buildBondedToAtomList(topology)
#     for res in matched_res_list:
#         template, matches = forcefield._getResidueTemplateMatches(res, bondedToAtom)
#         if template in map_template_to_pmd:
#             map_omm_res_to_pmd[res] = map_template_to_pmd[template]
#         else:
#             res_top, res_pos = oeommutils.oemol_to_openmmTop(map_omm_to_oe[res])
#             try:
#                 res_omm_system = forcefield.createSystem(res_top, rigidWater=False, constraints=None)
#                 res_pmd = parmed.openmm.load_topology(res_top, res_omm_system, xyz=res_pos)
#             except:
#                 raise ValueError("Error in the recognised residue parametrization {}".format(res))
#
#             map_template_to_pmd[template] = res_pmd
#             map_omm_res_to_pmd[res] = res_pmd
#
#     # print(map_omm_res_to_pmd)
#
#     # UnMatched Residue Parametrization
#     for ures in unmatched_res_list:
#
#         oe_mol = map_omm_to_oe[ures]
#
#         # Charge the unrecognized excipient
#         if not oequacpac.OEAssignCharges(oe_mol, oequacpac.OEAM1BCCCharges(symmetrize=True)):
#             raise ValueError("It was not possible to charge the extract residue: {}".format(ures))
#
#         if other_ff in ['Smirnoff99Frosst', 'OpenFF_1.0.0', 'OpenFF_1.1.0']:
#             oe_mol = oeommutils.sanitizeOEMolecule(oe_mol)
#
#         pmd = ParamMolStructure(oe_mol, other_ff, prefix_name="MOL" + '_' + ures.name)
#         ures_pmd = pmd.parameterize()
#
#         map_omm_res_to_pmd[ures] = ures_pmd
#
#     # Component Parmed Structure
#     component_pmd = parmed.Structure()
#     for res in topology.residues():
#         component_pmd += map_omm_res_to_pmd[res]
#
#     if len(component_pmd.atoms) != component_copy.NumAtoms():
#         raise ValueError(
#             "Component OE molecule and Component Parmed structure number of atoms mismatch {} vs {}".format(
#                 len(component_pmd.atoms), component_copy.NumAtoms()))
#
#     # Set the positions
#     oe_comp_coord_dic = component_copy.GetCoords()
#     comp_coords = np.ndarray(shape=(component_copy.NumAtoms(), 3))
#     for at_idx in oe_comp_coord_dic:
#         comp_coords[at_idx] = oe_comp_coord_dic [at_idx]
#
#     component_pmd.coordinates = comp_coords
#
#     return component_pmd


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
