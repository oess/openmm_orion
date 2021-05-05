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
from simtk.openmm import app
import numpy as np
from oeommtools import utils as oeommutils
from oeommtools import data_utils
from pdbfixer import PDBFixer
from pkg_resources import resource_filename
from simtk import unit


def hydrate(system, opt):
    """
    This function solvates the system by using PDBFixer

    Parameters:
    -----------
    system: OEMol molecule
        The system to solvate
    opt: python dictionary
        The parameters used to solvate the system

    Return:
    -------
    oe_mol: OEMol
        The solvated system
    """
    def BoundingBox(molecule):
        """
        This function calculates the Bounding Box of the passed
        molecule

        molecule: OEMol

        return: bb (numpy array)
            the calculated bounding box is returned as numpy array:
            [(xmin,ymin,zmin), (xmax,ymax,zmax)]
        """
        coords = [v for k, v in molecule.GetCoords().items()]
        np_coords = np.array(coords)
        min_coord = np_coords.min(axis=0)
        max_coord = np_coords.max(axis=0)
        bb = np.array([min_coord, max_coord])
        return bb

    # Create a system copy
    sol_system = system.CreateCopy()

    # Calculate system BoundingBox (Angstrom units)
    BB = BoundingBox(sol_system)

    # Estimation of the box cube length in A
    box_edge = 2.0 * opt['solvent_padding'] + np.max(BB[1] - BB[0])

    # BB center
    xc = (BB[0][0]+BB[1][0])/2.
    yc = (BB[0][1]+BB[1][1])/2.
    zc = (BB[0][2]+BB[1][2])/2.

    delta = np.array([box_edge/2., box_edge/2., box_edge/2.]) - np.array([xc, yc, zc])

    sys_coord_dic = {k: (v+delta) for k, v in sol_system.GetCoords().items()}

    sol_system.SetCoords(sys_coord_dic)

    # Load a fake system to initialize PDBfixer
    filename = resource_filename('pdbfixer', 'tests/data/test.pdb')
    fixer = PDBFixer(filename=filename)

    # Convert between OE and OpenMM topology
    omm_top, omm_pos = oeommutils.oemol_to_openmmTop(sol_system)

    chain_names = []

    for chain in omm_top.chains():
        chain_names.append(chain.id)

    # Set the correct topology to the fake system
    fixer.topology = omm_top
    fixer.positions = omm_pos

    # Solvate the system
    fixer.addSolvent(padding=unit.Quantity(opt['solvent_padding'], unit.angstroms),
                     ionicStrength=unit.Quantity(opt['salt_concentration'], unit.millimolar))

    # The OpenMM topology produced by the solvation fixer has missing bond
    # orders and aromaticity. The following section is creating a new openmm
    # topology made of just water molecules and ions. The new topology is then
    # converted in an OEMol and added to the passed molecule to produce the
    # solvated system

    wat_ion_top = app.Topology()

    # Atom dictionary between the the PDBfixer topology and the water_ion topology
    fixer_atom_to_wat_ion_atom = {}

    for chain in fixer.topology.chains():
        if chain.id not in chain_names:
            n_chain = wat_ion_top.addChain(chain.id)
            for res in chain.residues():
                n_res = wat_ion_top.addResidue(res.name, n_chain)
                for at in res.atoms():
                    n_at = wat_ion_top.addAtom(at.name, at.element, n_res)
                    fixer_atom_to_wat_ion_atom[at] = n_at

    for bond in fixer.topology.bonds():
        at0 = bond[0]
        at1 = bond[1]
        try:
            wat_ion_top.addBond(fixer_atom_to_wat_ion_atom[at0],
                                fixer_atom_to_wat_ion_atom[at1], type="Single", order=1)
        except:
            pass

    wat_ion_pos = fixer.positions[len(omm_pos):]

    oe_wat_ions = oeommutils.openmmTop_to_oemol(wat_ion_top, wat_ion_pos)

    oechem.OEAddMols(sol_system, oe_wat_ions)

    # Setting the box vectors
    omm_box_vectors = fixer.topology.getPeriodicBoxVectors()
    box_vectors = data_utils.encodePyObj(omm_box_vectors)
    sol_system.SetData(oechem.OEGetTag('box_vectors'), box_vectors)

    return sol_system


def order_check(mol, fname):
    """
    TO REMOVE
    This function is used to debug
    """
    import logging
    logger = logging.getLogger('Testing')
    hdlr = logging.FileHandler(fname)
    formatter = logging.Formatter('%(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)

    # hv = oechem.OEHierView(mol, oechem.OEAssumption_BondedResidue + oechem.OEAssumption_ResPerceived)

    hv = oechem.OEHierView(mol)

    for chain in hv.GetChains():
        logger.info('{}'.format(chain.GetChainID()))
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                logger.info('\t{} {}'.format(hres.GetOEResidue().GetName(), hres.GetOEResidue().GetResidueNumber()))
                for oe_at in hres.GetAtoms():
                    logger.info('\t\t{} {}'.format(oe_at.GetName(), oe_at.GetIdx()))

    return


def clash_detection(core_mol, check_mol, cutoff=2.0):

    def info(info_list):
        info_str = "Atom        Atom        Dist(A)  ratio\n"
        info_str += 38 * '-' + "\n"
        for tup in info_list:
            info_str += "{:<11s} {:<11s} {:1.5f}  {:1.3f}\n".format(tup[0].GetName() + " " + str(tup[0].GetIdx()),
                                                                    tup[1].GetName() + " " + str(tup[1].GetIdx()),
                                                                    tup[3],
                                                                    tup[2])

        return info_str

    # Create the Nearest neighbours
    nn = oechem.OENearestNbrs(check_mol, cutoff)

    severe_risk_list = list()
    moderate_risk_list = list()
    low_risk_list = list()

    # Check neighbours setting the atom bit mask
    for nbrs in nn.GetNbrs(core_mol):
        atom_bgn = nbrs.GetBgn()
        atom_end = nbrs.GetEnd()
        dist = nbrs.GetDist()

        if atom_bgn.GetAtomicNum() == 1:
            continue
        else:
            atom_bgn_vdwr = oechem.OEGetBondiVdWRadius(atom_bgn.GetAtomicNum())
        if atom_end.GetAtomicNum() == 1:
            continue
        else:
            atom_end_vdwr = oechem.OEGetBondiVdWRadius(atom_end.GetAtomicNum())

        # Polar Hydrogen VdW radii settings
        # if atom_bgn.IsPolarHydrogen():
        #     atom_bgn_vdwr = 0.95
        # else:
        #     atom_bgn_vdwr = oechem.OEGetBondiVdWRadius(atom_bgn.GetAtomicNum())
        # if atom_end.IsPolarHydrogen():
        #     atom_end_vdwr = 0.95
        # else:
        #     atom_end_vdwr = oechem.OEGetBondiVdWRadius(atom_end.GetAtomicNum())

        alpha = dist / (atom_bgn_vdwr + atom_end_vdwr)

        if alpha <= 0.3:
            severe_risk_list.append((atom_bgn, atom_end, alpha, dist))
        elif 0.3 < alpha <= 0.5:
            moderate_risk_list.append((atom_bgn, atom_end, alpha, dist))
        else:
            low_risk_list.append((atom_bgn, atom_end, alpha, dist))

    severe_risk_list.sort(key=lambda x: x[2])
    moderate_risk_list.sort(key=lambda x: x[2])
    low_risk_list.sort(key=lambda x: x[2])

    severe_info = info(severe_risk_list)
    moderate_info = info(moderate_risk_list)
    low_risk_info = info(low_risk_list)

    return [severe_risk_list, severe_info], [moderate_risk_list, moderate_info], [low_risk_list, low_risk_info]
