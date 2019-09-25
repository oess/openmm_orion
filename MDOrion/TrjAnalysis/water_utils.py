from openeye import (oechem,
                     oegrid,
                     oespicoli)


def dist2(coord1, coord2):
    distsq = (coord1[0] - coord2[0])**2 + \
             (coord1[1] - coord2[1])**2 + \
             (coord1[2] - coord2[2])**2

    return distsq


def nmax_waters(protein, ligand, cutoff):

    # Grid Spacing in A
    spacing = 0.5

    complex = oechem.OEMol(protein)

    oechem.OEAddMols(complex, ligand)

    surf = oespicoli.OESurface()
    oespicoli.OEMakeMolecularSurface(surf, complex, spacing)

    # oespicoli.OEWriteSurface("test_surf.oesrf", surf)

    center = oechem.OEFloatArray(3)
    extents = oechem.OEFloatArray(3)

    oechem.OEGetCenterAndExtents(ligand, center, extents)

    extents = oechem.OEFloatArray([max(extents) * 2, max(extents) * 2, max(extents) * 2])

    grid_reference = oegrid.OEScalarGrid()
    oegrid.OEMakeGridFromCenterAndExtents(grid_reference, center, extents, spacing)

    grid_spicoli = oegrid.OEScalarGrid()
    oespicoli.OEMakeBitGridFromSurface(grid_spicoli, surf)

    for iz in range(grid_reference.GetZDim()):
        for iy in range(grid_reference.GetYDim()):
            for ix in range(grid_reference.GetXDim()):
                x = grid_reference.GetX(ix)
                y = grid_reference.GetY(iy)
                z = grid_reference.GetZ(iz)

                ix_spicoli = grid_spicoli.GetXIdx(x)
                iy_spicoli = grid_spicoli.GetYIdx(y)
                iz_spicoli = grid_spicoli.GetZIdx(z)

                value = grid_spicoli.GetValue(ix_spicoli,
                                              iy_spicoli,
                                              iz_spicoli)

                grid_reference.SetValue(ix,
                                        iy,
                                        iz,
                                        value)
    # Invert Grid
    for iz in range(grid_reference.GetZDim()):
        for iy in range(grid_reference.GetYDim()):
            for ix in range(grid_reference.GetXDim()):
                # print("ix = {} iy = {} iz = {} value = {}".format(ix, iy, iz , grid.GetValue(ix, iy, iz)))
                if grid_reference.GetValue(ix, iy, iz) == 0.0:
                    grid_reference.SetValue(ix, iy, iz, 1.0)
                else:
                    grid_reference.SetValue(ix, iy, iz, 0.0)

    ligand_coords = ligand.GetCoords()

    for iz in range(grid_reference.GetZDim()):
        for iy in range(grid_reference.GetYDim()):
            for ix in range(grid_reference.GetXDim()):
                if grid_reference.GetValue(ix, iy, iz) == 1.0:

                    x = grid_reference.GetX(ix)
                    y = grid_reference.GetY(iy)
                    z = grid_reference.GetZ(iz)

                    min_sq = cutoff * cutoff

                    for coord in ligand_coords.values():
                        distsq = dist2(coord, (x, y, z))
                        if distsq < min_sq:
                            min_sq = distsq
                            break
                    if min_sq == cutoff * cutoff:
                        grid_reference.SetValue(ix, iy, iz, 0.0)

    protein_coords = protein.GetCoords()

    for iz in range(grid_reference.GetZDim()):
        for iy in range(grid_reference.GetYDim()):
            for ix in range(grid_reference.GetXDim()):
                if grid_reference.GetValue(ix, iy, iz) == 1.0:

                    x = grid_reference.GetX(ix)
                    y = grid_reference.GetY(iy)
                    z = grid_reference.GetZ(iz)

                    min_sq = cutoff * cutoff

                    for coord in protein_coords.values():
                        distsq = dist2(coord, (x, y, z))
                        if distsq < min_sq:
                            min_sq = distsq
                            break
                    if min_sq == cutoff * cutoff:
                        grid_reference.SetValue(ix, iy, iz, 0.0)

    # oegrid.OEWriteGrid("protein_ligand.grd", grid_reference)

    count = 0

    for iz in range(grid_reference.GetZDim()):
        for iy in range(grid_reference.GetYDim()):
            for ix in range(grid_reference.GetXDim()):
                if grid_reference.GetValue(ix, iy, iz) == 1.0:
                    count += 1

    # Calculate Volume from count in Angstrom
    vcount = (spacing ** 3) * count

    # Number of water molecule in vcount volume
    nwaters = int(0.034 * vcount)

    return nwaters


# def select_nmax_waters()