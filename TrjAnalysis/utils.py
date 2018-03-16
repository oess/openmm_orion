import numpy as np


def BoundingBox(molecule, save_pdb=False, scale_factor=1.0):
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

    min_coord_sf = np.array([min_coord[0] * (1.0 + scale_factor) + max_coord[0] * (1.0-scale_factor),
                             min_coord[1] * (1.0 + scale_factor) + max_coord[1] * (1.0 - scale_factor),
                             min_coord[2] * (1.0 + scale_factor) + max_coord[2] * (1.0 - scale_factor)])*0.5

    max_coord_sf = np.array([min_coord[0] * (1.0 - scale_factor) + max_coord[0] * (1.0 + scale_factor),
                             min_coord[1] * (1.0 - scale_factor) + max_coord[1] * (1.0 + scale_factor),
                             min_coord[2] * (1.0 - scale_factor) + max_coord[2] * (1.0 + scale_factor)])*0.5

    bb = np.array([min_coord_sf, max_coord_sf])

    if save_pdb:
        cube = [
                [min_coord_sf[0], min_coord_sf[1], min_coord_sf[2]],  # 0
                [max_coord_sf[0]-0.5, min_coord_sf[1], min_coord_sf[2]],  # 1
                [min_coord_sf[0], max_coord_sf[1]-0.5, min_coord_sf[2]],  # 2
                [max_coord_sf[0]-0.5, max_coord_sf[1]-0.5, min_coord_sf[2]],  # 3
                [min_coord_sf[0], min_coord_sf[1], max_coord_sf[2]-0.5],  # 4
                [max_coord_sf[0]-0.5, min_coord_sf[1], max_coord_sf[2]-0.5],  # 5
                [min_coord_sf[0], max_coord_sf[1]-0.5, max_coord_sf[2]-0.5],  # 6
                [max_coord_sf[0]-0.5, max_coord_sf[1]-0.5, max_coord_sf[2]-0.5]   # 7
                ]

        f = open("bb.pdb", 'w')
        count = 0
        at = 'H'
        bf = 1.0
        for coord in cube:
            ln = "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s  " % (count % 100000, at,
                                                                                            'CUB', 'A',
                                                                                            1,
                                                                                            coord[0],
                                                                                            coord[1],
                                                                                            coord[2],
                                                                                            bf,
                                                                                            'H')
            f.write(ln + '\n')
            count += 1

        center = "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s  " % (8 % 100000, 'C',
                                                                                            'CUB', 'A',
                                                                                            1,
                                                                                            0.5*(min_coord_sf[0]+max_coord_sf[0]),
                                                                                            0.5*(min_coord_sf[1]+max_coord_sf[1]),
                                                                                            0.5*(min_coord_sf[2]+max_coord_sf[2]),
                                                                                            bf,
                                                                                            'C')
        f.write(center + '\n')

        connect = ['CONECT    0    1    2    4',
                   'CONECT    1    0    3    5',
                   'CONECT    2    0    3    6',
                   'CONECT    3    1    2    7',
                   'CONECT    3    1    2    7',
                   'CONECT    4    5    6    0',
                   'CONECT    5    4    7    1',
                   'CONECT    6    4    7    2',
                   'CONECT    7    5    6    3']

        for ct in connect:
            f.write(ct + '\n')

        f.close()

    return bb
