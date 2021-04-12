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

import numpy as np

from openeye import (oechem,
                     oedepict,
                     oegrapheme,
                     oegrid,
                     oespicoli,
                     oedocking)
import mdtraj as md

from datarecord import OEField

import os

import contextlib

import glob

from tempfile import TemporaryDirectory

import oetrajanalysis.Clustering_utils as clusutl

from MDOrion.TrjAnalysis.water_utils import nmax_waters


def extract_aligned_prot_lig_wat_traj(md_components, flask, trj_fn, opt, nmax=30, water_cutoff=15.0):
    """
    Extracts the aligned protein trajectory and aligned ligand trajectory and aligned
    Water trajectory from a MD trajectory of a larger system that includes other
    components (eg water).
    The passed in setup mol must have the topology that matches the trajectory, and its xyz
    coordinates are the reference for the alignment. The alignment is done on the
    alpha carbons (atom name CA) of the active site residues within cutoff
    from the ligand. Once the alignment is done, the protein and ligand trajectories
    are each placed into a separate OEMol, one conformer per trajectory frame.
    Water trajectory is selecting the nmax waters from the ligand and protein CA
    within the cutoff distance for each trajectory snapshot

    Inputs:
        md_components: MDComponents object
            The md components carrying the setup starting flask.

        flask: OEMol
            The system flask

        trj_fn: String
            The filename of the hdf5-format MD trajectory or Gromacs .trr file format
        water_cutoff: Float
            The cutoff distance between the PL binding site and the waters in angstroms
        nmax: Integer
            max number of waters to select
    Outputs:
        multi_conf_protein: A multi conformer OEMol for the protein, one conformer per frame.
        multi_conf_ligand: A multi conformer OEMol for the ligand, one conformer per frame.
        multi_conf_water: A multi conformer OEMol for the waters, one conformer per frame.
    """

    # Extract protein, ligand, water and excipients from the flask
    # protein, ligand, water, excipients = oeommutils.split(flask, ligand_res_name="LIG")

    set_up_flask, map_dic = md_components.create_flask
    protein = md_components.get_protein
    ligand = md_components.get_ligand

    check_nmax = nmax_waters(protein, ligand, water_cutoff)

    if check_nmax < nmax:
        opt['Logger'].warn("The selected number of max waters cannot fit around the protein binding site: {} vs {}".
                           format(nmax, check_nmax))

    void, traj_ext = os.path.splitext(trj_fn)

    traj_dir = os.path.dirname(trj_fn)

    if traj_ext == '.h5':
        trj = md.load_hdf5(trj_fn)

    elif traj_ext == '.trr':
        pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
        trj = md.load_trr(trj_fn, top=pdb_fn)
        trj = trj[1:]
    else:
        raise ValueError("Trajectory file format {} not recognized in the trajectory {}".format(traj_ext, trj_fn))

    # System topology
    top_trj = trj.topology

    # Ligand indexes
    # lig_idx = top_trj.select("resname LIG")
    lig_idx = map_dic['ligand']

    # Protein indexes
    # prot_idx = top_trj.select("protein")

    # It is safer to use OE toolkits than mdtraj which is missing the protein caps
    prot_idx = map_dic['protein']

    # for at in protein.GetAtoms():
    #     prot_idx.append(at.GetIdx())

    # Water oxygen indexes
    water_O_idx = top_trj.select("water and element O")

    # Protein carbon alpha indexes
    prot_ca_idx = top_trj.select("backbone and element C")

    # Cutoff for the selection of the binding site atoms in A
    cutoff_bs = 5.0

    # Carbon alpha binding site indexes
    ca_bs_idx = md.compute_neighbors(trj[0], cutoff_bs / 10.0, lig_idx, haystack_indices=prot_ca_idx, periodic=True)[0]

    # Carbon alpha binding site and ligand indexes
    ca_bs_lig_idx = np.concatenate((ca_bs_idx, lig_idx))

    # Image the protein-ligand trajectory so the complex does not jump across box boundaries
    protlig = trj[0].atom_slice(np.concatenate((prot_idx, lig_idx)))
    protligAtoms = [atom for atom in protlig.topology.atoms]

    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stderr(devnull):
            trjImaged = trj.image_molecules(inplace=False, anchor_molecules=[protligAtoms], make_whole=True)

    # trjImaged = trj.image_molecules(inplace=False, anchor_molecules=[protligAtoms], make_whole=True)

    count = 0
    water_max_frames = []

    # TODO DEBUG
    # trjImaged = trjImaged[:10]

    for frame in trjImaged:
        # print(count, flush=True)

        # Water oxygen binding site indexes
        water_O_bs_idx = md.compute_neighbors(frame,
                                              water_cutoff / 10.0,
                                              ca_bs_lig_idx,
                                              haystack_indices=water_O_idx,
                                              periodic=True)

        # Pair combination water indexes times ligand indexes
        wat_lig_pairs = np.array(np.meshgrid(water_O_bs_idx, lig_idx)).T.reshape(-1, 2)

        # Distances between the waters and the ligand in nm
        wat_lig_distances = md.compute_distances(frame, wat_lig_pairs, periodic=True, opt=True)

        # Reshape the wat_lig_distances
        ns = np.reshape(wat_lig_distances, (len(water_O_bs_idx[0]), len(lig_idx)))

        # Min distances in nm between the oxygen waters and the ligand
        min_wat_O_lig_distances = np.min(ns, axis=1)

        # Pair combination water indexes times protein binding site carbon alpha indexes
        wat_ca_bs_pairs = np.array(np.meshgrid(water_O_bs_idx, ca_bs_idx)).T.reshape(-1, 2)

        # Distances between the waters and the protein binding site carbon alpha in nm
        wat_ca_bs_distances = md.compute_distances(frame, wat_ca_bs_pairs, periodic=True, opt=True)

        # Reshape the wat_ca_bs_distances
        ns = np.reshape(wat_ca_bs_distances, (len(water_O_bs_idx[0]), len(ca_bs_idx)))

        # Min distances in nm between the oxygen waters and the protein binding site carbon alpha
        min_wat_O_ca_bs_distances = np.min(ns, axis=1)

        metrics = min_wat_O_lig_distances + min_wat_O_ca_bs_distances

        metric_distances = dict()

        for wat_idx, m in zip(water_O_bs_idx[0], metrics):
            metric_distances[int(wat_idx)] = m

        water_list_sorted_max = sorted(metric_distances.items(), key=lambda x: x[1])[:nmax]

        if len(water_list_sorted_max) != nmax:
            raise ValueError("The ordered water list has the wrong size {} vs expected {} for the frame {}".
                             format(len(water_list_sorted_max), nmax, count))

        water_max_frames.append(water_list_sorted_max)

        # print(min_wat_O_ca_bs_distances)
        # print(pairs[:len(lig_idx), :])
        # for p,d in zip(wat_ca_bs_pairs, wat_ca_bs_distances[0]):
        #     print(p,d)

        count += 1

    # Put the reference mol xyz into the 1-frame topologyTraj to use as a reference in the fit
    setup_mol_array_coords = oechem.OEDoubleArray(3 * set_up_flask.GetMaxAtomIdx())
    set_up_flask.GetCoords(setup_mol_array_coords)

    setup_mol_xyzArr = np.array(setup_mol_array_coords)
    setup_mol_xyzArr.shape = (-1, 3)

    trj_reference = trjImaged[0]
    # convert from angstroms to nanometers
    trj_reference.xyz[0] = setup_mol_xyzArr / 10.0

    # Fitting
    trjImaged.superpose(trj_reference, 0, ca_bs_idx)

    # Delete Original Trajectory to save memory
    del trj

    # Molecule copies
    ligand_reference = oechem.OEMol(ligand)
    protein_reference = oechem.OEMol(protein)

    count = 0

    # Create the multi conformer protein, ligand and water molecules
    for frame in trjImaged.xyz:
        # print("Trj Image loop", count, flush=True)

        # Extract coordinates in A
        xyz = frame * 10

        # Set flask Coordinates as the current frame for the water extraction
        flask.SetCoords(xyz.flatten())
        water_list_sorted_max = water_max_frames[count]

        # print(water_list_sorted_max)

        # TODO The following solution to extract the waters do not
        #  keep the water order

        # Mark the close water atoms and extract them
        bv = oechem.OEBitVector(flask.GetMaxAtomIdx())
        bv.ClearBits()

        water_idx = []

        for pair in water_list_sorted_max:

            ow = flask.GetAtom(oechem.OEHasAtomIdx(pair[0]))

            # Select the whole water molecule
            for atw in oechem.OEGetResidueAtoms(ow):
                bv.SetBitOn(atw.GetIdx())
                water_idx.append(atw.GetIdx())

        pred_vec = oechem.OEAtomIdxSelected(bv)
        water_nmax_reference = oechem.OEMol()
        oechem.OESubsetMol(water_nmax_reference, flask, pred_vec)

        # TODO The following solution to extract the waters
        #  keep the water order but is it seems extremely inefficient

        # water_list = []
        # for pair in water_list_sorted_max:
        #     bv = oechem.OEBitVector(3)
        #     water_idx = []
        #     ow = flask.GetAtom(oechem.OEHasAtomIdx(pair[0]))
        #
        #     # Select the whole water molecule
        #     for atw in oechem.OEGetResidueAtoms(ow):
        #         bv.SetBitOn(atw.GetIdx())
        #         water_idx.append(atw.GetIdx())
        #
        #     pred_vec = oechem.OEAtomIdxSelected(bv)
        #     water = oechem.OEMol()
        #     oechem.OESubsetMol(water, flask, pred_vec)
        #
        #     water_list.append(water)
        #
        #
        # # print(len(water_list))
        #
        # water_nmax_reference = oechem.OEMol()

        # for w in water_list:
        #     oechem.OEAddMols(water_nmax_reference, w)

        # ligand and protein conf coordinates
        lig_xyz_list = [10 * frame[idx] for idx in lig_idx]
        lig_confxyz = oechem.OEFloatArray(np.array(lig_xyz_list).ravel())

        prot_xyz_list = [10 * frame[idx] for idx in prot_idx]
        prot_confxyz = oechem.OEFloatArray(np.array(prot_xyz_list).ravel())

        # Initialize the protein, ligand and water molecule topologies
        if count == 0:

            multi_conf_water = oechem.OEMol(water_nmax_reference)

            if multi_conf_water.NumAtoms() % 3 != 0:
                raise ValueError("Number of Water atoms is not multiple of 3")

            # Clean ResNumber and Chain on the multi conf water molecule
            # oechem.OEPerceiveResidues(multi_conf_water, oechem.OEPreserveResInfo_All)
            multi_conf_water.SetTitle("Water_" + str(nmax))

            res_num = 0
            i = 0
            for at in multi_conf_water.GetAtoms():

                res = oechem.OEAtomGetResidue(at)
                res.SetSerialNumber(i)
                res.SetName("HOH")
                res.SetChainID("Z")
                if i % 3 == 0:
                    res_num += 1
                res.SetResidueNumber(res_num)
                i += 1

            ligand_reference.SetCoords(lig_confxyz)
            protein_reference.SetCoords(prot_confxyz)
            multi_conf_ligand = oechem.OEMol(ligand_reference)
            multi_conf_protein = oechem.OEMol(protein_reference)

        # Attach the conformers on the multi conformer protein, ligand and water molecules
        else:
            water_confxyz = oechem.OEFloatArray(water_nmax_reference.NumAtoms() * 3)
            water_nmax_reference.GetCoords(water_confxyz)

            multi_conf_water.NewConf(water_confxyz)
            multi_conf_ligand.NewConf(lig_confxyz)
            multi_conf_protein.NewConf(prot_confxyz)

        count += 1

    return multi_conf_protein, multi_conf_ligand, multi_conf_water


def RequestOEField(record, field, rType):
    if not record.has_value(OEField(field,rType)):
        # opt['Logger'].warn('Missing record field {}'.format( field))
        print('Missing record field {}'.format(field))
        raise ValueError('The record does not have field {}'.format(field))
    else:
        return record.get_value(OEField(field, rType))


def RequestOEFieldType(record, field):
    if not record.has_value(field):
        # opt['Logger'].warn('Missing record field {}'.format( field))
        print('Missing record field {}'.format(field.get_name()))
        raise ValueError('The record does not have field {}'.format(field.get_name()))
    else:
        return record.get_value(field)


def PoseInteractionsSVG(md_components, width=400, height=300):
    """Generate a OEGrapheme interaction plot for a protein-ligand complex.
    The input protein may have other non-protein components as well so
    the input protein is first split into components to isolate the protein
    only for the plot. This may have to be changed if other components need
    to be included in the plot.
    """
    # # perceive residue hierarchy of total system
    # if not oechem.OEHasResidues(proteinOrig):
    #     oechem.OEPerceiveResidues(proteinOrig, oechem.OEPreserveResInfo_All)
    #     print('Perceiving residues')

    # split the total system into components
    #protein, ligandPsuedo, water, other = oeommutils.split(proteinOrig)

    protein = md_components.get_protein
    ligand = md_components.get_ligand

    # make the OEHintInteractionContainer
    asite = oechem.OEInteractionHintContainer(protein, ligand)
    if not asite.IsValid():
        oechem.OEThrow.Fatal("Cannot initialize active site!")
    # do the perceiving
    oechem.OEPerceiveInteractionHints(asite)
    # set the depiction options
    opts = oegrapheme.OE2DActiveSiteDisplayOptions(width, height)
    opts.SetRenderInteractiveLegend(True)
    magnifyresidue = 1.0
    opts.SetSVGMagnifyResidueInHover(magnifyresidue)
    # make the depiction
    oegrapheme.OEPrepareActiveSiteDepiction(asite)
    adisp = oegrapheme.OE2DActiveSiteDisplay(asite, opts)
    # make the image
    image = oedepict.OEImage(width, height)
    oegrapheme.OERenderActiveSite(image, adisp)
    # Add a legend
    #iconscale = 0.5
    #oedepict.OEAddInteractiveIcon(image, oedepict.OEIconLocation_TopRight, iconscale)
    svgBytes = oedepict.OEWriteImageToString("svg", image)

    svgString = svgBytes.decode("utf-8")

    return svgString


def GetBFactorColorGradient():
    colorg = oechem.OELinearColorGradient()
    colorg.AddStop(oechem.OEColorStop(0.0, oechem.OEDarkBlue))
    colorg.AddStop(oechem.OEColorStop(10.0, oechem.OELightBlue))
    colorg.AddStop(oechem.OEColorStop(25.0, oechem.OEYellowTint))
    colorg.AddStop(oechem.OEColorStop(50.0, oechem.OERed))
    colorg.AddStop(oechem.OEColorStop(100.0, oechem.OEDarkRose))
    return colorg


def ligand_to_svg_stmd(ligand, ligand_name):

    class ColorLigandAtomByBFactor(oegrapheme.OEAtomGlyphBase):
        def __init__(self, colorg):
            oegrapheme.OEAtomGlyphBase.__init__(self)
            self.colorg = colorg

        def RenderGlyph(self, disp, atom):
            adisp = disp.GetAtomDisplay(atom)
            if adisp is None or not adisp.IsVisible():
                return False

            if not oechem.OEHasResidue(atom):
                return False

            res = oechem.OEAtomGetResidue(atom)
            bfactor = res.GetBFactor()
            color = self.colorg.GetColorAt(bfactor)

            pen = oedepict.OEPen(color, color, oedepict.OEFill_On, 1.0)
            radius = disp.GetScale() / 3.0

            layer = disp.GetLayer(oedepict.OELayerPosition_Below)
            circlestyle = oegrapheme.OECircleStyle_Default
            oegrapheme.OEDrawCircle(layer, circlestyle, adisp.GetCoords(), radius, pen)
            return True

        def CreateCopy(self):
            return ColorLigandAtomByBFactor(self.colorg).__disown__()

    with TemporaryDirectory() as output_directory:

        lig_copy = oechem.OEMol(ligand)

        if len(ligand_name) < 15:
            lig_copy.SetTitle(ligand_name)
        else:
            lig_copy.SetTitle(ligand_name[0:13] + '...')

        img_fn = os.path.join(output_directory, "img.svg")

        oegrapheme.OEPrepareDepictionFrom3D(lig_copy)

        b_factor_color_grad = GetBFactorColorGradient()
        color_bfactor = ColorLigandAtomByBFactor(b_factor_color_grad)

        width, height = 150, 150
        opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
        opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)
        disp = oedepict.OE2DMolDisplay(lig_copy, opts)

        oegrapheme.OEAddGlyph(disp, color_bfactor, oechem.OEIsTrueAtom())

        oedepict.OERenderMolecule(img_fn, disp)

        svg_lines = ""
        marker = False
        with open(img_fn, 'r') as file:
            for line in file:
                if marker:
                    svg_lines += line

                if line.startswith("<svg"):
                    marker = True
                    svg_lines += line
                    svg_lines += """<title>{}</title>\n""".format(ligand_name)

                if line.startswith("</svg>"):
                    marker = False

    return svg_lines


def HighlightStyleMolecule(mol):
    hiliteColorer = oechem.OEMolStyleColorer(oechem.OEAtomColorScheme_Element)
    hiliteCarbonColor = oechem.OEColor( 245, 210,150 )
    hiliteColorer.AddColor(6, hiliteCarbonColor )
    hiliteColorer.AddColor(15, oechem.OEMagenta)
    #
    hiliteConfStyle = oechem.OE3DMolStyle()
    hiliteConfStyle.SetAtomStyle(oechem.OEAtomStyle_Stick)
    hiliteConfStyle.SetHydrogenVisibility(oechem.OEHydrogenVisibility_Polar)
    hiliteConfStyle.SetAtomColorer(hiliteColorer)
    #
    oechem.OEClearStyle( mol)
    oechem.OESetStyle( mol, hiliteConfStyle)
    return


def SetProteinLigandVizStyle(protein, ligand, carbonRGB=(180, 180, 180)):
    # set the carbon color (and phosphorus to magenta)
    carbonColor = oechem.OEColor(carbonRGB[0], carbonRGB[1], carbonRGB[2])
    acolorer = oechem.OEMolStyleColorer(oechem.OEAtomColorScheme_Element)
    acolorer.AddColor(15, oechem.OEMagenta)
    acolorer.AddColor(6, carbonColor)
    # make the protein style object
    protein_style = oechem.OE3DMolStyle()
    protein_style.SetHydrogenVisibility(oechem.OEHydrogenVisibility_Polar)
    protein_style.SetAtomStyle(oechem.OEAtomStyle_Hidden)
    protein_style.SetProteinStyle(oechem.OEProteinStyle_Ribbons)
    protein_style.SetProteinColorer(oechem.OEMolStyleColorer(oechem.OEProteinColorScheme_AtomColor))
    protein_style.SetAtomColorer(acolorer)
    # make the active site style object
    asite_style = oechem.OE3DMolStyle()
    asite_style.SetAtomStyle(oechem.OEAtomStyle_Wireframe)
    # make the ligand style object
    ligand_style = oechem.OE3DMolStyle()
    ligand_style.SetAtomStyle(oechem.OEAtomStyle_Stick)
    ligand_style.SetHydrogenVisibility(oechem.OEHydrogenVisibility_Polar)
    ligand_style.SetAtomColorer(acolorer)
    # now color the protein and the ligand
    oechem.OEClearStyle( protein)
    oechem.OESetStyle( protein, protein_style)
    oechem.OEClearStyle( ligand)
    oechem.OESetStyle( ligand, ligand_style)
    # display binding site residues
    asitePred = oechem.OEAtomMatchResidue( protein, ligand, 5)
    for atom in protein.GetAtoms( asitePred):
        oechem.OEClearStyle(atom)
        oechem.OESetStyle(atom, asite_style)
    return


def StyleTrajProteinLigandClusters( protein, ligand):
    if protein.NumConfs() != ligand.NumConfs():
        print('Cannot style; protein and ligand must have same number of conformers')
        return False
    confRGB = clusutl.ColorblindRGBMarkerColors( protein.NumConfs())
    # print( confRGB)
    SetProteinLigandVizStyle( protein, ligand, confRGB[0])
    for pconf, lconf, colorRGB in zip(protein.GetConfs(), ligand.GetConfs(), confRGB):
        # print( pconf.GetTitle(), lconf.GetTitle(), colorRGB)
        SetProteinLigandVizStyle( pconf, lconf, colorRGB)
    return True


class ClusterStyler:

    def __init__(self, num_clusters=0):
        self._confRGB = clusutl.ColorblindRGBMarkerColors(num_clusters)

    def apply_style(self, du, cluster_id):
        # TODO: Consider special style on waters and cofactors too
        SetProteinLigandVizStyle(du.GetImpl().GetProtein(), du.GetImpl().GetLigand(), self._confRGB[cluster_id])
        return True
