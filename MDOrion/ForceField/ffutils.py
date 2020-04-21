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


from openeye import (oechem,
                     oequacpac)

from oeommtools import utils as oeommutils

import numpy as np

import tempfile

import parmed

import openmoltools

from openmoltools.openeye import *

from simtk.openmm.app import AmberInpcrdFile, AmberPrmtopFile
from simtk.openmm import app

from openforcefield.typing.engines.smirnoff import ForceField

from pkg_resources import resource_filename

from openforcefield.topology import Topology, Molecule

from oeommtools.utils import oemol_to_openmmTop

from MDOrion.LigPrep.ff_utils import assignELF10charges

from MDOrion.ForceField.ff_library import ff_library


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

    def __init__(self, molecule, forcefield, prefix_name='ligand', delete_out_files=True):
        if not forcefield in list(ff_library.ligandff.values()):
            raise RuntimeError('The selected ligand force field is not '
                               'supported {}. Available {}'.format(forcefield, list(ff_library.ligandff.keys())))
        else:
            self.molecule = molecule
            self.forcefield = str(forcefield).strip()
            self.structure = None
            self.prefix_name = prefix_name
            self.delete_out_files = delete_out_files

    def checkCharges(self, molecule):
        # Check that molecule is charged.
        is_charged = False
        for atom in molecule.GetAtoms():
            if atom.GetPartialCharge() != 0.0:
                is_charged = True
        if not is_charged:
            raise Exception('Molecule {} has no charges; input molecules must be charged.' .format(molecule.GetTitle()))

    def getSmirnoffStructure(self, molecule=None):

        if not molecule:
            molecule = self.molecule
        try:
            self.checkCharges(molecule)
        except:
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")

            molecule = assignELF10charges(molecule)

        if self.forcefield == ff_library.ligandff['Smirnoff99Frosst']:

            fffn = resource_filename('openforcefield', os.path.join('data', 'test_forcefields/' + self.forcefield))

            if not os.path.exists(fffn):
                raise ValueError(
                    "Sorry! {} does not exist. If you just added it, you'll have to re-install".format(fffn))

            with open(fffn) as ffxml:
                ff = ForceField(ffxml, allow_cosmetic_attributes=True)

        elif self.forcefield in [ff_library.ligandff['OpenFF_1.0.0'],  ff_library.ligandff['OpenFF_1.1.0']]:
            ff = ForceField(self.forcefield, allow_cosmetic_attributes=True)

        else:
            raise ValueError("Force Field not Supported: {}".format(self.forcefield))

        mol_off = Molecule.from_openeye(molecule, allow_undefined_stereo=True)
        topology = Topology.from_molecules([mol_off])

        omm_sys = ff.create_openmm_system(topology, charge_from_molecules=[mol_off])

        # omm_top = generateTopologyFromOEMol(molecule)
        # positions = mol_off.conformers[0]

        omm_top, positions = oemol_to_openmmTop(molecule)

        pmd_structure = parmed.openmm.load_topology(omm_top, omm_sys, xyz=positions)

        return pmd_structure

    def getGaffStructure(self, molecule=None, forcefield=None):
        if not molecule:
            molecule = self.molecule
        if not forcefield:
            forcefield = self.forcefield

        if not molecule:
            molecule = self.molecule
        try:
            self.checkCharges(molecule)
        except:
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")

            molecule = assignELF10charges(molecule)

        # Write out mol to a mol2 file to process via AmberTools
        mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
        mol2filename = mol2file.name
        with oechem.oemolostream(mol2filename) as ofs:
            oechem.OEWriteConstMolecule(ofs, molecule)

        # Run antechamber to type and parmchk for frcmod
        # requires openmoltools 0.7.5 or later, which should be conda-installable via omnia
        gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber(self.prefix_name, mol2filename,
                                                                                 gaff_version=forcefield.lower(),
                                                                                 charge_method=None)

        # Run tleap using specified forcefield
        prmtop, inpcrd = openmoltools.amber.run_tleap(self.prefix_name, gaff_mol2_filename,
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

        if self.forcefield in [ff_library.ligandff['OpenFF_1.0.0'], ff_library.ligandff['OpenFF_1.1.0'], ff_library.ligandff['Smirnoff99Frosst']]:
            structure = self.getSmirnoffStructure()

        elif self.forcefield in ['GAFF', 'GAFF2']:
            structure = self.getGaffStructure()

        self.structure = structure
        return self.structure


def parametrize_component(component, component_ff):

    component_copy = oechem.OEMol(component)

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(component_copy)

    # Try to apply the selected FF on the component
    forcefield = app.ForceField(component_ff)

    # List of the unrecognized component
    unmatched_res_list = forcefield.getUnmatchedResidues(topology)

    if not unmatched_res_list:
        omm_components = forcefield.createSystem(topology, rigidWater=False, constraints=None)
        components_pmd = parmed.openmm.load_topology(topology, omm_components, xyz=positions)
        return components_pmd, None

    numparts, partlist = oechem.OEDetermineComponents(component_copy)
    pred = oechem.OEPartPredAtom(partlist)

    part_mols = []
    for i in range(1, numparts + 1):
        pred.SelectPart(i)
        partmol = oechem.OEMol()
        oechem.OESubsetMol(partmol, component_copy, pred)
        part_mols.append(partmol)

    map_omm_to_oe = {omm_res: oe_mol for omm_res, oe_mol in zip(topology.residues(), part_mols)}

    matched_res_list = []
    for res in topology.residues():
        if res in unmatched_res_list:
            continue
        else:
            matched_res_list.append(res)

    if not matched_res_list:
        return None, component_copy

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
                raise ValueError("Error in the recognised other residue parametrization {}".format(res))

            map_template_to_pmd[template] = res_pmd
            map_omm_res_to_pmd[res] = res_pmd

    # print(map_omm_res_to_pmd)

    # Matched OE Mol
    matched_oe_mols = oechem.OEMol()

    # Component Parmed Structure
    matched_pmd = parmed.Structure()

    for res in matched_res_list:
        matched_pmd += map_omm_res_to_pmd[res]
        oe_part = map_omm_to_oe[res]
        oechem.OEAddMols(matched_oe_mols, oe_part)

    if len(matched_pmd.atoms) != matched_oe_mols.NumAtoms():
        raise ValueError(
            "The OE Component molecule and the corresponding Parmed structure have "
            "number of atoms mismatch {} vs {}".format(
                matched_oe_mols.NumAtoms(), len(matched_pmd.atoms)))

    oe_comp_coord_dic = matched_oe_mols.GetCoords()
    comp_coords = np.ndarray(shape=(matched_oe_mols.NumAtoms(), 3))

    for at_idx in oe_comp_coord_dic:
        comp_coords[at_idx] = oe_comp_coord_dic[at_idx]

    matched_pmd.coordinates = comp_coords

    un_matched_oe_mols = oechem.OEMol()
    for un_res in unmatched_res_list:
        oechem.OEAddMols(un_matched_oe_mols, map_omm_to_oe[un_res])

    if un_matched_oe_mols.NumAtoms() == 0:
        return matched_pmd, None
    else:
        un_matched_oe_mols.SetTitle("Unmatched_Cofactors")
        return matched_pmd, un_matched_oe_mols


def parametrize_unknown_component(component, other_ff):

    component_copy = oechem.OEMol(component)

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(component_copy)

    numparts, partlist = oechem.OEDetermineComponents(component_copy)
    pred = oechem.OEPartPredAtom(partlist)

    part_mols = []
    for i in range(1, numparts + 1):
        pred.SelectPart(i)
        partmol = oechem.OEMol()
        oechem.OESubsetMol(partmol, component_copy, pred)
        part_mols.append(partmol)

    map_omm_to_oe = {omm_res: oe_mol for omm_res, oe_mol in zip(topology.residues(), part_mols)}

    map_omm_res_name_to_pmd = dict()
    for res in topology.residues():

        res_name = res.name

        if res_name not in map_omm_res_name_to_pmd:

            oe_mol = map_omm_to_oe[res]

            # Charge the unrecognized component
            if not oequacpac.OEAssignCharges(oe_mol, oequacpac.OEAM1BCCCharges(symmetrize=True)):
                raise ValueError("Is was not possible to charge the extract residue: {}".format(res))

            if other_ff in ['Smirnoff99Frosst', 'OpenFF_1.0.0', 'OpenFF_1.1.0']:
                oe_mol = oeommutils.sanitizeOEMolecule(oe_mol)

            pmd = ParamMolStructure(oe_mol, other_ff, prefix_name="MOL" + '_' + res_name)
            res_pmd = pmd.parameterize()

            map_omm_res_name_to_pmd[res_name] = res_pmd

    # Whole Parmed Structure
    comp_pmd = parmed.Structure()
    for res in topology.residues():
        res_name = res.name
        comp_pmd += map_omm_res_name_to_pmd[res_name]

    if len(comp_pmd.atoms) != component_copy.NumAtoms():
        raise ValueError(
            "The OE Component molecule and the corresponding Parmed structure have "
            "number of atoms mismatch {} vs {}".format(
                component_copy.NumAtoms(), len(comp_pmd.atoms)))

    # Set the positions
    oe_comp_coord_dic = component_copy.GetCoords()
    comp_coords = np.ndarray(shape=(component_copy.NumAtoms(), 3))
    for at_idx in oe_comp_coord_dic:
        comp_coords[at_idx] = oe_comp_coord_dic[at_idx]

    comp_pmd.coordinates = comp_coords

    return comp_pmd
