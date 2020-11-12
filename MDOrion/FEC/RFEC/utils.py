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


from pyparsing import (Word,
                       printables)

import os

import parmed

from tempfile import TemporaryDirectory

from oemdtoolbox.FEC.RBFEC.chimera import Chimera

import mdtraj as md

from openeye import oechem

import numpy as np

from simtk import unit

from oeommtools.utils import openmmTop_to_oemol

from MDOrion.FEC.RFEC.gmx_templates import gromacs_nes

import MDOrion.FEC.RFEC.pmx as pmx

import subprocess

from subprocess import STDOUT, PIPE, Popen, DEVNULL

import tarfile

from simtk.openmm import app

import glob


def edge_map_grammar(word):

    map_string = Word(printables) + ">>" + Word(printables)

    expr = map_string.parseString(word)

    # expr = word.rstrip().split()
    # print(expr)

    return expr


def parmed_find_ligand(pmd, lig_res_name="LIG"):
    pmd_split = pmd.split()

    # print(pmd_split)

    for idx in range(0, len(pmd_split)):
        pmd_struc = pmd_split[idx][0]

        for res in pmd_struc.residues:
            if res.name == lig_res_name:
                return pmd_struc, idx

    return None, None


def gmx_chimera_topology_injection(pmd_flask, pmd_chimera_start, pmd_chimera_final):

    with TemporaryDirectory() as outputdir:

        # outputdir = "./"

        flask_top_fn = os.path.join(outputdir, "flask.top")
        gmx_chimera_fn_prefix = os.path.join(outputdir, "gmx_chimera")

        top_gmx = parmed.gromacs.GromacsTopologyFile.from_structure(pmd_flask)
        top_gmx.defaults = parmed.gromacs.gromacstop._Defaults(fudgeLJ=0.5, fudgeQQ=0.8333, gen_pairs='yes')
        top_gmx.write(flask_top_fn)

        with open(flask_top_fn, 'r') as f:
            flask_gmx_topology_lines = f.readlines()

        gmx_chimera_atomtype, gmx_chimera_moltype = Chimera.to_gromacs(pmd_chimera_start, pmd_chimera_final,
                                                                       filename=gmx_chimera_fn_prefix,
                                                                       mixed_mode=True,
                                                                       itp=True,
                                                                       perturb_mass=True)

    # Merge the Atom Type Section between the Flask and the Chimera
    capture = False
    end_idx = 0

    for idx, ln in enumerate(flask_gmx_topology_lines):
        if '[ atomtypes ]' in ln:
            capture = True
            continue
        if ln.startswith('[') or ln.startswith("#"):
            capture = False
        if capture:
            end_idx = idx
            continue

    flask_gmx_topology_lines.insert(end_idx-1, "\n".join(gmx_chimera_atomtype.split("\n")[2:]))

    # Change molecule section to add the chimera gmx molecule type
    capture = False
    start_idx = 0
    end_idx = 0

    for idx, ln in enumerate(flask_gmx_topology_lines):
        if '[ moleculetype ]' in ln and 'LIG' in flask_gmx_topology_lines[idx + 2]:
            start_idx = idx
            capture = True
            continue
        elif '[ moleculetype ]' in ln:
            capture = False
        if capture:
            end_idx = idx
            continue

    del flask_gmx_topology_lines[start_idx:end_idx]

    flask_gmx_topology_lines.insert(start_idx, gmx_chimera_moltype)

    # Update the molecule section
    capture = False
    for idx, ln in enumerate(flask_gmx_topology_lines):
        if '[ molecules ]' in ln:
            capture = True
            continue
        if capture and ln.startswith("LIG"):
            flask_gmx_topology_lines[idx] = "CMR                  1\n"
            capture = False

    gmx_out = "".join(flask_gmx_topology_lines)

    return gmx_out


def gmx_chimera_coordinate_injection(pmd_chimera, mdrecord, tot_frames, query_mol, chimera, chimera_graph_dic):

    md_components = mdrecord.get_md_components
    set_up_flask, map_dic = md_components.create_flask
    ligand = md_components.get_ligand

    lig_idx = map_dic['ligand']

    trj_fn = mdrecord.get_stage_trajectory(stg_name='last')
    trj = md.load(trj_fn)

    # trj.save("trj.pdb")

    if len(trj) > tot_frames:
        stride = int(len(trj)/tot_frames)
        trj = trj[0:tot_frames*stride:stride]

    map_excess_final_new_idxs = chimera_graph_dic['map_excess_final_new_idxs']
    map_chimera_to_final_idxs = chimera_graph_dic['map_chimera_to_final_idxs']
    morph = chimera_graph_dic['direction']

    if morph == "A_to_B":
        pmd_initial = chimera.pmdA
    else:
        pmd_initial = chimera.pmdB

    pmd_flask = mdrecord.get_parmed(sync_stage_name="last")

    new_pmd_structure = parmed.Structure()
    before_lig_at = []
    after_lig_at = []
    switch = False
    for res in pmd_flask.residues:
        for at in res.atoms:
            if res.name == 'LIG':
                switch = True
                continue
            if not switch:
                before_lig_at.append(at.idx)
            else:
                after_lig_at.append(at.idx)

    ligand_reference = oechem.OEMol(ligand)

    gmx_gro_str_list = list()

    for count, tr in enumerate(trj):
        print(">>>", count)
        frame_xyz = tr.xyz[0] * 10
        bv = tr.unitcell_vectors[0] * 10 * unit.angstrom

        lig_xyz_list = [frame_xyz[idx] for idx in lig_idx]

        lig_confxyz = oechem.OEFloatArray(np.array(lig_xyz_list).ravel())
        ligand_reference.SetCoords(lig_confxyz)

        best_conf = chimera._alignment(ligand_reference, query_mol, morph)

        best_conf_coords = best_conf.GetCoords()

        pmd_initial.coordinates = np.array(lig_xyz_list).reshape(ligand_reference.NumAtoms(), 3)

        pmd_chimera_initial_coords = {}
        for at_chimera in pmd_chimera.atoms:
            if at_chimera.idx not in map_excess_final_new_idxs.values():
                pmd_chimera_initial_coords[at_chimera.idx] = pmd_initial.coordinates[at_chimera.idx]
            else:
                # print(">>>>>", at_chimera.idx, map_chimera_to_final_idxs[at_chimera.idx])
                coord = np.array(best_conf_coords[map_chimera_to_final_idxs[at_chimera.idx]])
                pmd_chimera_initial_coords[at_chimera.idx] = coord

        sorted_coords = np.array([p[1] for p in sorted(pmd_chimera_initial_coords.items())])
        pmd_chimera.coordinates = sorted_coords.reshape(len(pmd_chimera.atoms), 3)

        pmd_flask.coordinates = frame_xyz
        pmd_flask.box_vectors = bv

        if count == 0:
            if before_lig_at:
                new_pmd_structure += pmd_flask[before_lig_at]

            new_pmd_structure += pmd_chimera

            if after_lig_at:
                new_pmd_structure += pmd_flask[after_lig_at]

            new_pmd_structure.box_vectors = pmd_flask.box_vectors

            # TODO DEBUG ONLY MULTI CHIMERA
            # oe_chimera = openmmTop_to_oemol(pmd_chimera.topology, pmd_chimera.positions)

        else:

            if before_lig_at:
                max_idx = max(before_lig_at)
            else:
                max_idx = 0

            if max_idx > 0:
                pmd_flask_before_lig_coords = pmd_flask.coordinates[0:max_idx+1, :]
                pmd_flask_after_lig_coords = pmd_flask.coordinates[max_idx + ligand.NumAtoms() + 1:, :]
                new_coords = np.concatenate((pmd_flask_before_lig_coords,
                                             pmd_chimera.coordinates,
                                             pmd_flask_after_lig_coords))
            else:
                pmd_flask_after_lig_coords = pmd_flask.coordinates[ligand.NumAtoms():, :]
                new_coords = np.concatenate((pmd_chimera.coordinates, pmd_flask_after_lig_coords))

            if len(new_pmd_structure.atoms) != new_coords.shape[0]:
                raise ValueError("Coordinate atom mismatch size {} vs {}".format(len(new_pmd_structure.atoms),
                                                                                 new_coords.shape[0]))
            new_pmd_structure.box_vectors = pmd_flask.box_vectors
            new_pmd_structure.coordinates = new_coords

        with TemporaryDirectory() as outputdir:

            # outputdir = "./"

            flask_gro_fn = os.path.join(outputdir, "flask.gro")
            new_pmd_structure.save(flask_gro_fn, overwrite=True)

            with open(flask_gro_fn, 'r') as f:
                flask_gmx_gro_str = f.read()

            gmx_gro_str_list.append(flask_gmx_gro_str)

        # TODO DEBUG ONLY MULTI CHIMERA
        # pos = new_pmd_structure.coordinates
        # chimera_oe_pos = oechem.OEFloatArray(pos.ravel())
        # oe_chimera.NewConf(chimera_oe_pos)

    # TODO DEBUG ONLY MULTI CHIMERA
    # with oechem.oemolostream("chimera_multi.pdb") as ofs:
    #     oechem.OEWriteConstMolecule(ofs, oe_chimera)

    return gmx_gro_str_list


def gmx_nes_run(gmx_gro, gmx_top, opt):

    out_dir = opt['out_directory']

    gmx_gro_fn = os.path.join(out_dir, "gmx_gro.gro")
    gmx_top_fn = os.path.join(out_dir, "gmx_top.top")
    gmx_ne_mdp_fn = os.path.join(out_dir, "gmx_ne.mdp")
    gmx_ne_tpr_fn = os.path.join(out_dir, "gmx_ne.tpr")
    gmx_deffnm_out = os.path.join(out_dir, "gmx_run" + '_' + str(opt['frame_count']))
    gmx_trj_fn = gmx_deffnm_out + '.trr'

    with open(gmx_gro_fn, 'w') as f:
        f.write(gmx_gro)

    with open(gmx_top_fn, 'w') as f:
        f.write(gmx_top)

    stepLen = 2.0 * unit.femtoseconds
    nsteps = int(round(opt['time'] / (stepLen.in_units_of(unit.nanoseconds) / unit.nanoseconds)))
    if opt['enable_switching']:
        dlambda = 1.0/nsteps
    else:
        dlambda = 0

    min_box = opt['min_box']

    # Cutoff in A
    cutoff = 10

    # in A
    theshold = (min_box / 2.0) * 0.85

    if cutoff < theshold:
        cutoff_distance = cutoff * unit.angstroms
    else:
        cutoff_distance = theshold * unit.angstroms

        opt['Logger'].warn("Cutoff Distance too large for the box size. Set the cutoff distance "
                           "to {:.2f} A".format(cutoff_distance.value_in_unit(unit.angstrom)))

    rvdw_switch = cutoff_distance - 1.0 * unit.angstrom

    pressure = opt['pressure'] * unit.atmosphere

    gmx_fe_template = gromacs_nes.format(nsteps=nsteps,
                                         temperature=opt['temperature'],
                                         pressure=pressure.value_in_unit(unit.bar),
                                         cutoff=cutoff_distance.value_in_unit(unit.nanometer),
                                         rvdwswitch=rvdw_switch.value_in_unit(unit.nanometer),
                                         dlambda=dlambda)

    with open(gmx_ne_mdp_fn, 'w') as f:
        f.write(gmx_fe_template)

    if opt['verbose']:

        # Assemble the Gromacs system to run
        subprocess.check_call(['gmx',
                               'grompp',
                               '-f', gmx_ne_mdp_fn,
                               '-c', gmx_gro_fn,
                               '-p', gmx_top_fn,
                               '-o', gmx_ne_tpr_fn,
                               '-maxwarn', '4'
                               ])

        # Run Gromacs
        subprocess.check_call(['gmx',
                               'mdrun',
                               '-v',
                               '-s', gmx_ne_tpr_fn,
                               '-deffnm', gmx_deffnm_out
                               ])

        # Convert the trajectory in .xtc
        # p = subprocess.Popen(['gmx',
        #                       'trjconv',
        #                       '-f', gmx_trj_fn,
        #                       '-s', gmx_ne_tpr_fn,
        #                       '-o', gmx_deffnm_out + '.xtc',
        #                       '-pbc', b'whole'],
        #                      stdin=subprocess.PIPE)
        #
        # # Select the entire System
        # p.communicate(b'0')

    else:
        # Assemble the Gromacs system to run
        subprocess.check_call(['gmx',
                               'grompp',
                               '-f', gmx_ne_mdp_fn,
                               '-c', gmx_gro_fn,
                               '-p', gmx_top_fn,
                               '-o', gmx_ne_tpr_fn,
                               '-maxwarn', '4'
                               ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)

        # Run Gromacs
        subprocess.check_call(['gmx',
                               'mdrun',
                               '-v',
                               '-s', gmx_ne_tpr_fn,
                               '-deffnm', gmx_deffnm_out
                               ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)

        # Convert the trajectory in .xtc
        # p = subprocess.Popen(['gmx',
        #                       'trjconv',
        #                       '-f', gmx_trj_fn,
        #                       '-s', gmx_ne_tpr_fn,
        #                       '-o', gmx_deffnm_out + '.xtc',
        #                       '-pbc', b'whole'],
        #                      stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)
        #
        # # Select the entire System
        # p.communicate(b'0')

    opt['log_fn'] = gmx_deffnm_out + '.log'

    tar_fn = opt['trj_fn']

    with tarfile.open(tar_fn, mode='w:gz') as archive:
        archive.add(opt['out_directory'], arcname='.')

    import shutil
    shutil.copy(tar_fn, "/Users/gcalabro/Desktop")

    # Generate the final Parmed structure
    gro = app.GromacsGroFile(gmx_deffnm_out+'.gro')
    top = app.GromacsTopFile(gmx_top_fn, unitCellDimensions=gro.getUnitCellDimensions())

    for chain in top.topology.chains():
        chain.id = 'A'

    omm_system = top.createSystem(rigidWater=False, constraints=None)
    final_pmd = parmed.openmm.load_topology(top.topology, omm_system, xyz=gro.positions)
    final_pmd.box_vectors = omm_system.getDefaultPeriodicBoxVectors()
    oe_flask = openmmTop_to_oemol(top.topology, gro.positions)

    opt['pmd'] = final_pmd
    opt['flask'] = oe_flask

    with open(gmx_deffnm_out+'.gro', 'r') as f:
        gro_str = f.read()

    opt['gro_str'] = gro_str

    if opt['enable_switching']:
        xvg_fn = glob.glob(os.path.join(out_dir, '*.xvg'))[0]

        w = pmx.parse_dgdl_files([xvg_fn], lambda0=0, invert_values=False)

        if w is None:
            raise ValueError("Work calculation failed")
        opt['gmx_work'] = w[0]
    else:
        opt['gmx_work'] = None

    return
