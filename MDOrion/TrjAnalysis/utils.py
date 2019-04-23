#############################################################################
# Copyright (C) 2019 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np

from openeye import (oechem,
                     oedepict,
                     oegrapheme)
import mdtraj as md

from datarecord import OEField

import os

import glob

from tempfile import TemporaryDirectory

from MDOrion.Standards import MDEngines

def GetCardinalOrderOfProteinResNums(mol):
    # make map of protein res nums to the residue cardinal order index
    resmap = {}
    currRes = -10000
    currIdx = -1
    for atom in mol.GetAtoms(oechem.OEIsBackboneAtom()):
        thisRes = oechem.OEAtomGetResidue(atom)
        resnum = thisRes.GetResidueNumber()
        if resnum != currRes:
            currIdx += 1
            currRes = resnum
            resmap[currRes] = currIdx
    return resmap, currIdx


def ExtractProtLigActsiteResNums(mol, fromLigCutoff=5.0):
    '''Extracts the protein and ligand from a single OEMol containing a protein-ligand
    complex plus other components. A list of protein residues within a cutoff distance
    from the ligand is also returned.
    Inputs:
        mol: The OEMol containing protein, ligand, and all other components
        fromLigCutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
    Returns:
        protein: An OEMol containing only the protein
        ligand: An OEMol containing only the ligand
        actSiteResNums: A list of integers, one per residue number for protein residues
            with any atom within fromLigCutoff distance of the ligand.'''
    # perceive residue hierarchy of total system
    if not oechem.OEHasResidues(mol):
        oechem.OEPerceiveResidues(mol, oechem.OEPreserveResInfo_All)
    # split the total system into components
    ligand = oechem.OEMol()
    protein = oechem.OEMol()
    water = oechem.OEMol()
    other = oechem.OEMol()
    #sopts = oechem.OESplitMolComplexOptions('MOL')
    sopts = oechem.OESplitMolComplexOptions()
    oechem.OESplitMolComplex(ligand, protein, water, other, mol, sopts)
    # use residue-based distance cutoff (in angstroms) from ligand to define an active site residue
    nn = oechem.OENearestNbrs(protein, fromLigCutoff)
    actSiteResNums = set()
    for nbrs in nn.GetNbrs(ligand):
        residue = oechem.OEAtomGetResidue(nbrs.GetBgn())
        actSiteResNums.add(residue.GetResidueNumber())
    return protein, ligand, actSiteResNums


def ExtractAlignedProtLigTraj(mol, traj_filename, fromLigCutoff=5.0, skip=0):
    '''Extracts the aligned protein trajectory and aligned ligand trajectory from
    a MD trajectory of a larger system that includes other components (eg water).
    The passed in OEMol must have the topology that matches the trajectory, and its xyz
    coordinates are the reference for the alignment. The alignment is done on the
    alpha carbons (atom name CA) of the active site residues within fromLigCutoff
    from the ligand. Once the alignment is done, the protein and ligand trajectories
    are each placed into a separate OEMol, one conformer per trajectory frame.
    Inputs:
        mol: An OEMol giving the topology for the trajectory and the reference xyz
            coordinates for the alignment.
        traj_filename: The filename of the hdf5-format MD trajectory or Gromacs .xtc file format
        fromLigCutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
        skip: number of frames to skip at the beginning of the trajectory.
    Outputs:
        protTraj: A multiconformer OEMol for the protein, one conformer per frame.
        ligTraj: A multiconformer OEMol for the ligand, one conformer per frame.'''

    void, traj_ext = os.path.splitext(traj_filename)

    traj_dir = os.path.dirname(traj_filename)

    # get the topology from 1st frame of the traj file
    if traj_ext == '.h5':
        topologyTraj = md.load_hdf5(traj_filename, frame=1)

    elif traj_ext == '.xtc':
        pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
        topologyTraj = md.load_xtc(traj_filename, top=pdb_fn, frame=1)
    else:
        raise ValueError("Trajectory file format {} not recognized in the trajecotry {}".format(traj_ext, traj_filename))

    # Put the reference mol xyz into the 1-frame topologyTraj to use as a reference in the fit
    molXyz = oechem.OEDoubleArray(3*mol.GetMaxAtomIdx())
    mol.GetCoords(molXyz)
    molXyzArr = np.array(molXyz)
    molXyzArr.shape = (-1, 3)

    # convert from angstroms to nanometers and slice out the protein-ligand complex
    topologyTraj.xyz[0] = molXyzArr/10.0

    # extract protein and ligand molecules from the larger multicomponent system
    # and identify residue numbers for residues within fromLigCutoff of the ligand.
    protein, ligand, actSiteResNums = ExtractProtLigActsiteResNums(mol, fromLigCutoff)
    protResMap, numProtRes = GetCardinalOrderOfProteinResNums(protein)
    actSiteResIdxs = set()
    for resnum in actSiteResNums:
        actSiteResIdxs.add( protResMap[resnum])

    # Extract protein atom indices: cannot trust mdtraj protein selection so
    # assume they are contiguous and starting the atom list and just get the same
    # number of atoms as in the OpenEye protein
    protOEIdx = np.array([atom.GetIdx() for atom in protein.GetAtoms()])

    # Extract ligand atom indices
    #   Note: the ligand must have residue name 'MOL' or 'LIG' (bad, should change)
    ligIdx = topologyTraj.topology.select('resname == MOL or resname == LIG')
    protligIdx = np.append( protOEIdx, ligIdx)

    #print( 'numAtoms prot, lig, protlig:', len(protOEIdx), len(ligIdx), len(protligIdx))
    #protligIdx = topologyTraj.topology.select('protein or resname == MOL or resname == LIG')

    # Read the protein-ligand subsystem of the trajectory file
    if traj_ext == '.h5':  # OpenMM
        trj_initial = md.load_hdf5(traj_filename, atom_indices=protligIdx)
    elif traj_ext == '.xtc':  # Gromacs
        trj_initial = md.load_xtc(traj_filename, top=pdb_fn, atom_indices=protligIdx)

    if skip > 0 and len(trj_initial) > skip:
        trj = trj_initial[skip:]
    else:
        trj = trj_initial

    # Image the protein-ligand trajectory so the complex does not jump across box boundaries
    protlig = topologyTraj.atom_slice(protligIdx)
    protligAtoms = [atom for atom in protlig.topology.atoms]
    inplace = True
    trjImaged = trj.image_molecules(inplace, [protligAtoms])

    # Make a list of the atom indices of the carbon-alphas of the active site residues;
    # assume residue numbering matches the mol
    actSiteCA = [atom.index for atom in topologyTraj.topology.atoms
                 if ((atom.residue.resSeq in actSiteResIdxs) and (atom.name == 'CA'))]

    # Fit the protein-ligand trajectory to the active site carbon-alphas of the reference
    trjImaged.superpose(protlig, 0, actSiteCA)

    # Generate a multiconformer representation of the ligand trajectory
    ligTraj = oechem.OEMol(ligand)
    ligTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [10*frame[idx] for idx in ligIdx]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel())
        conf = ligTraj.NewConf(confxyz)

    # Generate a multiconformer representation of the protein trajectory
    strNumProteinAtomsToSelect = 'index ' + str(protOEIdx[0]) + ' to ' + str(protOEIdx[-1])
    protIdx = protlig.topology.select(strNumProteinAtomsToSelect)
    protTraj = oechem.OEMol(protein)
    protTraj.DeleteConfs()

    for frame in trjImaged.xyz:
        xyzList = [10*frame[idx] for idx in protIdx]
        confxyz = oechem.OEFloatArray(np.array(xyzList).ravel())
        conf = protTraj.NewConf(confxyz)

    return protTraj, ligTraj


def RequestOEField(record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))


def RequestOEFieldType( record, field):
    if not record.has_value(field):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field.get_name() ))
        raise ValueError('The record does not have field {}'.format( field.get_name() ))
    else:
        return record.get_value(field)


def ColorblindRGBMarkerColors( nColors=0):
    palette = [(0, 114, 178), (0, 158, 115), (213, 94, 0), (204, 121, 167),
        (240, 228, 66), (230, 159, 0), (86, 180, 233), (150, 150, 150)]
    if nColors<1:
        return palette
    elif nColors<9:
        return palette[:nColors]
    else:
        n = int(nColors/8)
        moreRGB = palette
        for i in range(n):
            moreRGB = moreRGB+palette
        return( moreRGB[:nColors])


def PoseInteractionsSVG(ligand, proteinOrig, width=400, height=300):
    """Generate a OEGrapheme interaction plot for a protein-ligand complex.
    The input protein may have other non-protein components as well so
    the input protein is first split into components to isolate the protein
    only for the plot. This may have to be changed if other components need
    to be included in the plot.
    """
    # perceive residue hierarchy of total system
    if not oechem.OEHasResidues(proteinOrig):
        oechem.OEPerceiveResidues(proteinOrig, oechem.OEPreserveResInfo_All)
        print('Perceiving residues')
    # split the total system into components
    ligandPsuedo = oechem.OEMol()
    protein = oechem.OEMol()
    water = oechem.OEMol()
    other = oechem.OEMol()
    #sopts = oechem.OESplitMolComplexOptions('MOL')
    sopts = oechem.OESplitMolComplexOptions()
    oechem.OESplitMolComplex(ligandPsuedo, protein, water, other, proteinOrig, sopts)
    #
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

    # TODO BUG TEMPORARY FIX FOR TOOLKIT VERSION 2018.10.1 JIRA CASE https://openeye.atlassian.net/browse/TOOLKS-624

    lines = svgString.splitlines(True)

    new_str = """.oedepict-visible {
     visibility: visible;
    }"""

    for idx in range(0, len(lines)):
        if lines[idx] == ' visibility: hidden;\n':
            lines.insert(idx + 2, new_str)
            break

    svg_out = " ".join(lines)

    # TODO END

    return svg_out


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

        colorg = oechem.OELinearColorGradient()
        colorg.AddStop(oechem.OEColorStop(0.0, oechem.OEDarkBlue))
        colorg.AddStop(oechem.OEColorStop(10.0, oechem.OELightBlue))
        colorg.AddStop(oechem.OEColorStop(25.0, oechem.OEYellowTint))
        colorg.AddStop(oechem.OEColorStop(50.0, oechem.OERed))
        colorg.AddStop(oechem.OEColorStop(100.0, oechem.OEDarkRose))

        color_bfactor = ColorLigandAtomByBFactor(colorg)

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
