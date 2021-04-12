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
                       alphanums,
                       pyparsing_common,
                       Literal,
                       Combine,
                       Optional,
                       printables)
import os

import parmed

from tempfile import TemporaryDirectory

from oemdtoolbox.FEC.RBFEC.chimera import Chimera

import mdtraj as md

from mdtraj.formats import TRRTrajectoryFile

from openeye import oechem

import numpy as np

from simtk import unit

from oeommtools.utils import openmmTop_to_oemol

import MDOrion.FEC.RFEC.pmx as pmx

from simtk.openmm import app

from simtk import openmm


from MDOrion.FEC.RFEC.pmx import BAR

import math

from openeye import oedepict

from plotly.subplots import make_subplots

import plotly

import plotly.graph_objects as go

from scipy.stats import norm

from PIL import Image

import h5py

from MDOrion.Standards.standards import Fields

from datarecord import Meta

import random

from MDOrion.FEC.RFEC import stats

from MDOrion.FEC.RFEC.freeenergyframework import wrangle

import scipy as sc

from MDOrion.FEC.RFEC.gmx_run import gmx_run

from orionclient.session import (in_orion,
                                 OrionSession,
                                 get_session)

from orionclient.types import (ShardCollection,
                               Shard)

from orionclient.helpers.collections import (try_hard_to_create_shard,
                                             try_hard_to_download_shard)

from MDOrion.Standards.standards import CollectionsNames

import glob

import subprocess


def edge_map_grammar(word):

    map_string = Word(printables) + ">>" + Word(printables)

    expr = map_string.parseString(word)

    # expr = word.rstrip().split()
    # print(expr)

    return expr


def rbfe_file_grammar(word):

    units = Combine(Literal("kcal/mol") | Literal("kJ/mol"))

    lig_name = Word(printables)

    DG = pyparsing_common.number

    dg_err = Optional(pyparsing_common.number)

    map_string = lig_name + DG + dg_err + units

    expr = map_string.parseString(word)

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


def fix_gromacs_water_names(pmd_flask):

    for r in pmd_flask.residues:

        if len(r.atoms) == 3:

            try:
                oxy, = (a for a in r.atoms if a.atomic_number == 8)

                hyd1, hyd2 = (a for a in r.atoms if a.atomic_number == 1)

                if oxy and hyd1 and hyd2:
                    oxy.name = 'O'
                    hyd1.name = 'H1'
                    hyd2.name = 'H2'
            except:
                continue

    return pmd_flask


def unique_atom_types(pmd_structure, ligand_res_name='LIG'):

    omm_system = pmd_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                            constraints=None,
                                            removeCMMotion=False,
                                            rigidWater=False)

    # # Copy the topology and positions
    topology = pmd_structure.topology
    positions = pmd_structure.positions

    def check_water(res):
        two_bonds = list(res.bonds())

        if len(two_bonds) == 2:

            waters = []

            for bond in two_bonds:

                elem0 = bond[0].element
                elem1 = bond[1].element

                if (elem0.atomic_number == 1 and elem1.atomic_number == 8) \
                        or (elem0.atomic_number == 8 and elem1.atomic_number == 1):
                    waters.append(True)

            if all(waters):
                return True

        else:
            return False

    atom_types_dic = {}
    count_id = 0

    for c in topology.chains():
        for r in c.residues():
            for a in r.atoms():
                if r.name + a.name in atom_types_dic:
                    a.id = atom_types_dic[r.name + a.name]
                else:
                    if check_water(r):
                        if a.element.atomic_number == 1:
                            a.id = 'HW'
                        else:
                            a.id = 'OW'
                        atom_types_dic[r.name + a.name] = a.id
                    elif r.name == ligand_res_name:
                        a.id = 'L' + str(count_id)
                    else:
                        a.id = 'O' + str(count_id)
                        atom_types_dic[r.name + a.name] = a.id

                    count_id += 1

    # Define a new parmed structure with the new unique atom types
    new_system_structure = parmed.openmm.load_topology(topology,
                                                       system=omm_system,
                                                       xyz=positions)
    # Re-set positions, velocities and box
    new_system_structure.positions = pmd_structure.positions
    new_system_structure.velocities = pmd_structure.velocities
    new_system_structure.box = pmd_structure.box

    return new_system_structure


def gmx_chimera_topology_injection(pmd_flask, pmd_chimera_start, pmd_chimera_final):

    pmd_flask = fix_gromacs_water_names(pmd_flask)

    pmd_flask = unique_atom_types(pmd_flask)

    with TemporaryDirectory() as outputdir:

        # TODO DEBUG
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
    # Clean Ligand Atom Type section
    capture = False
    idx_to_del = list()
    for idx, ln in enumerate(flask_gmx_topology_lines):
        if '[ atomtypes ]' in ln:
            capture = True
            continue
        if ln.startswith('[') or ln.startswith("#"):
            capture = False
        if capture:
            if ln.startswith('L'):
                idx_to_del.append(idx)

    flask_gmx_topology_lines_clean = list()
    for idx, ln in enumerate(flask_gmx_topology_lines):
        if idx not in idx_to_del:
            flask_gmx_topology_lines_clean.append(ln)

    # Merge the Atom Type Section between the Flask and the Chimera
    capture = False
    end_idx = 0

    for idx, ln in enumerate(flask_gmx_topology_lines_clean):
        if '[ atomtypes ]' in ln:
            capture = True
            continue
        if ln.startswith('[') or ln.startswith("#"):
            capture = False
        if capture:
            end_idx = idx
            continue

    flask_gmx_topology_lines_clean.insert(end_idx-1, "\n".join(gmx_chimera_atomtype.split("\n")[2:]))

    # Change molecule section to add the chimera gmx molecule type
    capture = False
    start_idx = 0
    end_idx = 0

    for idx, ln in enumerate(flask_gmx_topology_lines_clean):
        if '[ moleculetype ]' in ln and 'LIG' in flask_gmx_topology_lines_clean[idx + 2]:
            start_idx = idx
            capture = True
            continue
        elif '[ moleculetype ]' in ln:
            capture = False
        if capture:
            end_idx = idx
            continue

    del flask_gmx_topology_lines_clean[start_idx:end_idx]

    flask_gmx_topology_lines_clean.insert(start_idx, gmx_chimera_moltype)

    # Update the molecule section
    capture = False
    for idx, ln in enumerate(flask_gmx_topology_lines_clean):
        if '[ molecules ]' in ln:
            capture = True
            continue
        if capture and ln.startswith("LIG"):
            flask_gmx_topology_lines_clean[idx] = "CMR                  1\n"
            capture = False

    gmx_out = "".join(flask_gmx_topology_lines_clean)

    return gmx_out


def extract_trj_velocities(trj_fn,  indexes=None, trj_type="OpenMM"):

    if trj_type == "OpenMM":

        f = h5py.File(trj_fn, 'r')

        # velocities = f.get('velocities').value
        velocities = f['velocities'][()]

        f.close()
    elif trj_type == "Gromacs":
        trr = TRRTrajectoryFile(trj_fn, "r")
        v, f = trr._check_has_velocities_forces()
        if not v:
            raise ValueError("Gromacs Velocities missing")
        xyz, time, step, box, lam, velocities, force = trr._read(len(trr), indexes, get_velocities=True)

    else:
        raise ValueError("MD Engine not currently supported: {}".format(trj_type))

    return velocities


def gmx_chimera_coordinate_injection(pmd_chimera, mdrecord, tot_frames, query_mol, chimera, chimera_graph_dic):

    md_components = mdrecord.get_md_components
    set_up_flask, map_dic = md_components.create_flask
    ligand = md_components.get_ligand

    lig_idx = map_dic['ligand']

    trj_fn = mdrecord.get_stage_trajectory(stg_name='last')

    stage = mdrecord.get_stage_by_name('last')
    trj_field = stage.get_field(Fields.trajectory.get_name())
    trj_meta = trj_field.get_meta()
    md_engine = trj_meta.get_attribute(Meta.Annotation.Description)

    if md_engine == "OpenMM":
        trj = md.load(trj_fn)

        # Velocities array shape = (frame, atoms, 3) units nm/ps
        velocities = extract_trj_velocities(trj_fn, trj_type=md_engine)

    elif md_engine == "Gromacs":

        traj_dir = os.path.dirname(trj_fn)
        pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
        tpr_fn = glob.glob(os.path.join(traj_dir, '*.tpr'))[0]

        # To avoid issues with pbc conditions make whole the molecules
        trj_whole_fn = os.path.join(traj_dir, "trj_whole.trr")
        p = subprocess.Popen(['gmx',
                              'trjconv',
                              '-f', trj_fn,
                              '-s', tpr_fn,
                              '-o', trj_whole_fn,
                              '-pbc', b'whole'],
                             stdin=subprocess.PIPE)
        # Select the entire System
        p.communicate(b'0')

        trj = md.load_trr(trj_whole_fn, top=pdb_fn)

        # Skip first frame which is the loaded pdb topology
        trj = trj[1:]

        indexes = []
        for a in trj.topology.atoms:
            indexes.append(a.index)

        # Velocities array shape = (frame, atoms, 3) units nm/ps
        velocities = extract_trj_velocities(trj_whole_fn, indexes=indexes, trj_type=md_engine)

    else:
        raise ValueError("The selected md engine is not currently support".format(md_engine))

    if len(trj) > tot_frames:
        stride = int(len(trj)/tot_frames)
        trj = trj[0:tot_frames*stride:stride]

    map_excess_final_new_idxs = chimera_graph_dic['map_excess_final_new_idxs']
    map_chimera_to_final_idxs = chimera_graph_dic['map_chimera_to_final_idxs']
    morph = chimera_graph_dic['direction']

    if morph == "A_to_B":
        pmd_initial = chimera.pmdA
        pmd_final = chimera.pmdB
    else:
        pmd_initial = chimera.pmdB
        pmd_final = chimera.pmdA

    pmd_flask = mdrecord.get_parmed(sync_stage_name="last")

    # Fix water names
    pmd_flask = fix_gromacs_water_names(pmd_flask)

    # Generate Velocities from MB distributions to be assigned just to the
    # Dummy particles
    omm_sys_chimera = pmd_chimera.createSystem(nonbondedMethod=app.NoCutoff,
                                               constraints=None,
                                               removeCMMotion=False,
                                               rigidWater=False)
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    omm_chimera_context = openmm.Context(omm_sys_chimera, integrator)

    omm_chimera_context.setPositions(pmd_chimera.positions)
    omm_chimera_context.setVelocitiesToTemperature(300.0, random.randint(1, 100))

    omm_chimera_state = omm_chimera_context.getState(getVelocities=True)
    # Units nm/ps
    omm_chimera_state_vel = omm_chimera_state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)

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

    if before_lig_at:
        max_idx = max(before_lig_at)
    else:
        max_idx = 0

    ligand_reference = oechem.OEMol(ligand)

    gmx_gro_str_list = list()

    # TODO DEBUGGING
    for count, tr in enumerate(trj):
        # print(">>>", count)
        frame_xyz = tr.xyz[0] * 10
        bv = tr.unitcell_vectors[0] * 10 * unit.angstrom
        vel = velocities[count] * unit.nanometer / unit.picosecond

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

        # Optimize in place the pmd_chimera coordinates
        optmize_initial_atom_idxs = list(map_excess_final_new_idxs.values())
        chimera._optimize_chimera_coords(pmd_chimera, pmd_final, map_chimera_to_final_idxs, optmize_initial_atom_idxs)

        pmd_flask.coordinates = frame_xyz
        pmd_flask.box_vectors = bv
        pmd_flask.velocities = vel.in_units_of(unit.angstrom / unit.picosecond)

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
                raise ValueError("Atom Coordinates mismatch size {} vs {}".format(len(new_pmd_structure.atoms),
                                                                                  new_coords.shape[0]))
            new_pmd_structure.box_vectors = pmd_flask.box_vectors
            new_pmd_structure.coordinates = new_coords

        # Set Velocities
        cut_idx = len(pmd_chimera.atoms) - ligand.NumAtoms()
        dummy_velocities = omm_chimera_state_vel[0:cut_idx, :]
        if max_idx > 0:
            pmd_flask_before_lig_velocities = pmd_flask.velocities[0:max_idx + 1, :]
            pmd_flask_after_lig_velocities = pmd_flask.velocities[max_idx + ligand.NumAtoms() + 1:, :]
            lig_velocities = pmd_flask.velocities[max_idx + 1: max_idx + 1 + ligand.NumAtoms(), :]
            new_chimera_velocities = np.concatenate((lig_velocities, dummy_velocities))
            new_velocities = np.concatenate((pmd_flask_before_lig_velocities,
                                             new_chimera_velocities,
                                             pmd_flask_after_lig_velocities))

        else:
            lig_velocities = pmd_flask.velocities[0:ligand.NumAtoms(), :]
            new_chimera_velocities = np.concatenate((lig_velocities, dummy_velocities))
            pmd_flask_after_lig_velocities = pmd_flask.velocities[ligand.NumAtoms():, :]
            new_velocities = np.concatenate((new_chimera_velocities, pmd_flask_after_lig_velocities))

        new_pmd_structure.velocities = new_velocities

        if len(new_pmd_structure.velocities) != new_velocities.shape[0]:
            raise ValueError("Atom Velocities mismatch size {} vs {}".format(len(new_pmd_structure.velocities),
                                                                             new_velocities.shape[0]))

        with TemporaryDirectory() as outputdir:

            # outputdir = "./"

            flask_gro_fn = os.path.join(outputdir, "flask.gro")

            # TODO PARMED ERROR
            try:
                new_pmd_structure.save(flask_gro_fn, overwrite=True, combine='all')
                # new_pmd_structure.save(flask_gro_fn, overwrite=True)
            except:
                lig, lig_idx = parmed_find_ligand(new_pmd_structure, lig_res_name="CMR")

                split = new_pmd_structure.split()

                prot_max_idx = 0
                pmd_prot = parmed.Structure()
                for idx, p in enumerate(split):
                    if idx < lig_idx:
                        if len(p[1]) == 1:
                            prot_max_idx += len(p[0].atoms)
                            pmd_prot += p[0]
                        else:
                            raise ValueError("Parmed protein idexes Error")

                pmd_flask = pmd_prot + new_pmd_structure[prot_max_idx:]
                pmd_flask.box_vectors = new_pmd_structure.box_vectors
                pmd_flask.save(flask_gro_fn, overwrite=True)

            if not os.path.isfile(flask_gro_fn):
                raise ValueError("Gromacs Coordinate file Error")

            #  ###### DEBUGGING MISSING VELOCITIES WITH THE PDB FILE
            # flask_pdb_fn = os.path.join(outputdir, "flask.pdb")
            # new_pmd_structure.save(flask_pdb_fn, overwrite=True)
            # subprocess.check_call(['gmx',
            #                        'editconf',
            #                        '-f', flask_pdb_fn,
            #                        '-o', flask_gro_fn])
            #  ################### END DEBUGGING ###############

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

    gmx_deffnm_out = gmx_run(gmx_gro, gmx_top, opt)

    opt['log_fn'] = gmx_deffnm_out + '.log'

    tar_fn = opt['trj_fn']

    # TODO FOR DEBUGGING ONLY
    # import tarfile
    # with tarfile.open(tar_fn, mode='w:gz') as archive:
    #     archive.add(opt['out_directory'], arcname='.')

    # Generate the final Parmed structure
    gro = app.GromacsGroFile(gmx_deffnm_out+'.gro')
    top = app.GromacsTopFile(os.path.join(opt['out_directory'], "gmx_top.top"), unitCellDimensions=gro.getUnitCellDimensions())

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

    xvg_fn = gmx_deffnm_out + '.xvg'

    w = pmx.parse_dgdl_files([xvg_fn], lambda0=0, invert_values=False)

    if w is None:
        raise ValueError("Work calculation failed")

    opt['gmx_work'] = w[0]

    return


def nes_data_analysis(f_bound, r_bound, f_unbound, r_unbound):

    def init(f_dic, r_dic):

        frames_f = np.array([k for k, v in f_dic.items()])
        frames_r = np.array([k for k, v in r_dic.items()])
        work_f = np.array([v for k, v in f_dic.items()])
        work_r = np.array([v for k, v in r_dic.items()])

        return frames_f, frames_r, work_f, work_r

    results_analysis = dict()

    nboots = 100
    nblocks = 1

    # Extract data from the frame:work dictionary
    bound_frames_f, bound_frames_r, bound_work_f, bound_work_r = init(f_bound, r_bound)

    unbound_frames_f, unbound_frames_r, unbound_work_f, unbound_work_r = init(f_unbound, r_unbound)

    # BAR Analysis
    bar_bound = BAR(wf=bound_work_f, wr=bound_work_r, T=300.0, nboots=nboots, nblocks=nblocks)
    bar_unbound = BAR(wf=unbound_work_f, wr=unbound_work_r, T=300.0, nboots=nboots, nblocks=nblocks)

    binding_affinity = bar_bound.dg - bar_unbound.dg
    binding_affinity_err = math.sqrt(bar_bound.err_boot**2 + bar_unbound.err_boot**2)

    results_analysis['BAR'] = (binding_affinity, binding_affinity_err,
                               bar_bound.dg, bar_bound.err_boot,
                               bar_unbound.dg, bar_unbound.err_boot)
    #
    #     # Units are in kJ/mol
    #     print(">>>>>>>>>>>>>>", cgi.mf, cgi.mr, cgi.dg)
    #
    #     plot_work_dist(wf=res_ab, wr=res_ba, dG=cgi.dg)

    return results_analysis


def make_edge_depiction(ligandA, ligandB):

    ligandA_title = ligandA.GetTitle()
    ligandB_title = ligandB.GetTitle()

    for atom in ligandA.GetAtoms():
        atom.SetRxnRole(oechem.OERxnRole_Reactant)

    for atom in ligandB.GetAtoms():
        atom.SetRxnRole(oechem.OERxnRole_Product)

    reaction = oechem.OEMol()
    reaction.SetRxn(True)

    oechem.OEAddMols(reaction, ligandA)
    oechem.OEAddMols(reaction, ligandB)

    img_reaction_mol = oechem.OEMol(reaction)
    img_reaction_mol.SetTitle("")

    reaction.SetTitle(ligandA_title + " to " + ligandB_title)

    oedepict.OEPrepareDepiction(reaction)

    width, height = 500, 300
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    oe_disp = oedepict.OE2DMolDisplay(reaction, opts)

    with TemporaryDirectory() as outdir:

        fn = os.path.join(outdir, "depiction.svg")
        oedepict.OERenderMolecule(fn, oe_disp)

        with open(fn, 'r') as f:
            edge_depiction_string = f.read()

    oedepict.OEPrepareDepiction(img_reaction_mol)

    image = oedepict.OEImage(500, 300, oechem.OETransparentColor)
    clearbackground = True
    oedepict.OERenderMolecule(image, img_reaction_mol, not clearbackground)

    return edge_depiction_string, image


def plot_work_pdf(f_bound, r_bound, f_unbound, r_unbound, results, title, edge_depiction_string):

    def init(f_dic, r_dic):

        frames_f = np.array([k for k, v in f_dic.items()])
        frames_r = np.array([k for k, v in r_dic.items()])
        work_f = np.array([v for k, v in f_dic.items()])
        work_r = np.array([v for k, v in r_dic.items()])

        return frames_f, frames_r, work_f, work_r

    def make_plots(frames_f, frames_r, work_f, work_r, figure, row):

        # Forward color
        color_f = 'rgb(255,0,0)'

        # Reverse color
        color_r = 'rgb(0,0,255)'

        # Axis definition
        x_left = 'x'+str(row)
        x_right = 'x'+str(row+1)
        y_left = 'y'+str(row)
        y_right = 'y'+str(row+1)

        legend = False

        if row == 1:
            legend = True

        # Works vs Frames
        figure.add_trace(go.Scatter(x=frames_f, y=work_f,
                                    name='Forward Work',
                                    line=dict(color=color_f, width=2),
                                    mode='lines+markers',
                                    xaxis=x_left,
                                    yaxis=y_left,
                                    showlegend=False), row=row, col=1)

        figure.add_trace(go.Scatter(x=frames_r, y=work_r,
                                    name='Reverse Work',
                                    line=dict(color=color_r, width=2),
                                    mode='lines+markers',
                                    xaxis=x_right,
                                    yaxis=y_right,
                                    showlegend=False), row=row, col=1)

        # Histograms
        figure.add_trace(go.Histogram(y=work_f,
                                      histnorm='probability density',
                                      autobinx=True,
                                      opacity=0.75,
                                      marker_color=color_f,
                                      xaxis=x_right,
                                      yaxis=y_right,
                                      name="Forward Work",
                                      showlegend=legend,
                                      legendgroup="group"), row=row, col=2)

        figure.add_trace(go.Histogram(y=work_r,
                                      histnorm='probability density',
                                      autobinx=True,
                                      opacity=0.75,
                                      marker_color=color_r,
                                      xaxis=x_right,
                                      yaxis=y_right,
                                      name="Reverse Work",
                                      showlegend=legend,
                                      legendgroup="group", ), row=row, col=2)

        # Normal distributions
        mean_wf, sd_wf = norm.fit(work_f)
        mean_wr, sd_wr = norm.fit(work_r)

        works = np.concatenate((work_f, work_r))

        min_works = np.amin(works)
        max_works = np.amax(works)

        curve_wf = np.linspace(min_works, max_works, 500)
        pdf_wf = norm.pdf(curve_wf, loc=mean_wf, scale=sd_wf)

        curve_wr = np.linspace(min_works, max_works, 500)
        pdf_wr = norm.pdf(curve_wr, loc=mean_wr, scale=sd_wr)

        fig.add_trace(go.Scatter(x=pdf_wf, y=curve_wf,
                                 name='Forward Work',
                                 line=dict(color=color_f, width=2),
                                 mode='lines',
                                 xaxis=x_right,
                                 yaxis=y_right,
                                 showlegend=False), row=row, col=2)

        fig.add_trace(go.Scatter(x=pdf_wr, y=curve_wr,
                                 name='Forward Work',
                                 line=dict(color=color_r, width=2),
                                 mode='lines',
                                 xaxis=x_right,
                                 yaxis=y_right,
                                 showlegend=False), row=row, col=2)

        return figure

    # Extract data from the frame:work dictionary
    bound_frames_f, bound_frames_r, bound_work_f, bound_work_r = init(f_bound, r_bound)

    unbound_frames_f, unbound_frames_r, unbound_work_f, unbound_work_r = init(f_unbound, r_unbound)

    fig = make_subplots(rows=3, cols=2,
                        horizontal_spacing=0.01,
                        vertical_spacing=0.1,
                        column_widths=[0.7, 0.3],
                        shared_yaxes=True)
    # Bound Plot
    make_plots(bound_frames_f, bound_frames_r, bound_work_f, bound_work_r, fig, row=1)

    # Unbound Plot
    make_plots(unbound_frames_f, unbound_frames_r, unbound_work_f, unbound_work_r, fig, row=2)

    # Result table
    data = []
    methods = []
    ddg = []
    dgbound = []
    dgunbound = []

    for met, fecs in results.items():
        methods.append(met)
        ddg.append("{:.2f} \u00B1 {:.2f}".format(fecs[0], fecs[1]))
        dgbound.append("{:.2f} \u00B1 {:.2f}".format(fecs[2], fecs[3]))
        dgunbound.append("{:.2f} \u00B1 {:.2f}".format(fecs[4], fecs[5]))

    data.append(methods)
    data.append(ddg)
    data.append(dgbound)
    data.append(dgunbound)

    fig.add_trace(
        go.Table(
            header=dict(
                values=["Method", "\u0394\u0394G kJ/mol", "\u0394G Bound kJ/mol", "\u0394G Unbound kJ/mol"],
                font=dict(size=14),
                align="center",
                fill_color='paleturquoise',
            ),
            cells=dict(
                values=data,
                align="center"),

            domain=dict(x=[0, 1],
                        y=[0, 0.25])
        )

    )

    # Bound data Frames vs Work
    xaxis1 = dict(
        zeroline=True,
        showgrid=True,
        title="Frames"
    )

    yaxis1 = dict(
        zeroline=True,
        showgrid=True,
        title="Work Bound kJ/mol"
    )

    # Bound data PDF
    xaxis2 = dict(
        zeroline=False,
        showgrid=False,
        visible=True,
        title="PDF"
    )

    yaxis2 = dict(
        zeroline=False,
        showgrid=False,
        visible=False
    )

    # Bound data Frames vs Work
    xaxis3 = dict(
        zeroline=True,
        showgrid=True,
        title="Frames"
    )

    yaxis3 = dict(
        zeroline=True,
        showgrid=True,
        title="Work Unbound kJ/mol"
    )

    # Bound data PDF
    xaxis4 = dict(
        zeroline=False,
        showgrid=False,
        title="PDF"
    )

    yaxis4 = dict(
        zeroline=False,
        showgrid=False,
        visible=False
    )

    fig.update_layout(

        legend=dict(x=0.7, y=1.08),

        barmode='overlay',
        hovermode="closest",
        bargap=0,
        title="Non Equilibrium Switching Results: {}".format(title),

        # Bound data frames
        xaxis1=xaxis1,
        yaxis1=yaxis1,

        # Bound data PDF
        xaxis2=xaxis2,
        yaxis2=yaxis2,

        # Unbound data frames
        xaxis3=xaxis3,
        yaxis3=yaxis3,

        # Unbound data PDF
        xaxis4=xaxis4,
        yaxis4=yaxis4,

        # plot_bgcolor="rgba(0, 0, 0, 0)",
        # paper_bgcolor="rgba(0, 0, 0, 0)",
        autosize=False,
        width=900,
        height=900,
    )

    fig.update_yaxes(automargin=True)
    fig.update_xaxes(automargin=True)

    # fig.show(config={'scrollZoom': True})

    # fig.write_image("fig1.png")

    with TemporaryDirectory() as outdir:

        fn = os.path.join(outdir, "report.html")

        plotly.offline.plot(fig, filename=fn)

        with open(fn, 'r') as f:
            report_string_list = f.readlines()

        # Add the following lines for Orion CSS Plotly integration
        report_str = ""

        edge_string = """
        <style>
            .mutation {
                display: block;
                width: 400px;
                height: 250px;
        }

        .mutation svg {
            display: block;
            max-width: 100%;
            }
        </style>"""

        edge_string += """
        <main class="mutation">
            {svg}
        </main>\n
        """.format(svg=edge_depiction_string)

        for idx in range(0, len(report_string_list)):
            line = report_string_list[idx]
            if "<head>" in line:
                report_str += """<head><meta charset="utf-8" />   
                <style>
                    .modebar { display:unset; }
                    .js-plotly-plot { min-height: 900px; }
                    .js-plotly-plot { min-width: 900px; }
                </style></head>\n"""

            report_str += line

        report_str = edge_string + "\n" + report_str

    return report_str


def generate_plots_and_stats(exp_data_dic, predicted_data_dic, method='BAR', DDG_symmetrize=False, units='kcal/mol'):

    def calculate_statistics(exp, pred, plot_type='ddG'):

        try:
            mae = stats.bootstrap_statistic(exp, pred, statistic="MAE", plot_type=plot_type)
            rae = stats.bootstrap_statistic(exp, pred, statistic="RAE", plot_type=plot_type)
            rmse = stats.bootstrap_statistic(exp, pred, statistic="RMSE", plot_type=plot_type)
            rrmse = stats.bootstrap_statistic(exp, pred, statistic="RRMSE", plot_type=plot_type)
            r2 = stats.bootstrap_statistic(exp, pred, statistic="R2", plot_type=plot_type)
            rho = stats.bootstrap_statistic(exp, pred, statistic="RHO", plot_type=plot_type)
            ktau = stats.bootstrap_statistic(exp, pred, statistic="KTAU", plot_type=plot_type)

            # Generate statistic without bootstraping
            mae['data'] = stats.compute_statistic(exp, pred, "MAE")
            rae['data'] = stats.compute_statistic(exp, pred, "RAE")
            rmse['data'] = stats.compute_statistic(exp, pred, "RMSE")
            rrmse['data'] = stats.compute_statistic(exp, pred, "RRMSE")
            r2['data'] = stats.compute_statistic(exp, pred, "R2")
            rho['data'] = stats.compute_statistic(exp, pred, "RHO")
            ktau['data'] = stats.compute_statistic(exp, pred, "KTAU")

            statistics = dict()

            statistics['MAE'] = mae
            statistics['RAE'] = rae
            statistics['RMSE'] = rmse
            statistics['RRMSE'] = rrmse
            statistics['R2'] = r2
            statistics['RHO'] = rho
            statistics['KTAU'] = ktau

            return statistics

        except Exception as e:
            return None

    def sub_plot(x, err_x, y, err_y, figure, hover_text, col, range):

        slope, intercept, r_value, p_value, std_err = sc.stats.linregress(np.array(x),
                                                                          np.array(y))

        x_plt = np.linspace(range[0], range[1], 100, endpoint=True)
        line = slope * np.array(x_plt) + intercept

        x_plus_1kcal = x_plt + 4.184
        x_minus_1kcal = x_plt - 4.184

        alpha_color = 'rgba(127, 166, 238, 0.4)'

        # Theoretical line minus 1kcal
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=x_minus_1kcal,
            fill=None,
            hovertext="",
            hoverinfo="text",
            mode="lines",
            showlegend=False,
            line=dict(color=alpha_color, width=1)),
            row=1, col=col
        )

        # Theoretical line plus 1kcal
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=x_plus_1kcal,
            fill='tonexty',
            fillcolor=alpha_color,
            hovertext="",
            hoverinfo="text",
            mode="lines",
            showlegend=False,
            line=dict(color=alpha_color, width=1)),
            row=1, col=col
        )

        # Theoretical line
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=x_plt,
            hovertext="1:1",
            hoverinfo="text",
            mode="lines",
            showlegend=False,
            line=dict(color='black', width=2, dash='dash')),
            row=1, col=col
        )

        # Best fit line
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=line,
            hovertext="y = {:.2f} * x + {:.2f}".format(slope, intercept),
            hoverinfo="text",
            mode="lines",
            showlegend=False,
            line=dict(color='red', width=2)),
            row=1, col=col
        )

        if col == 1:
            hov_label = "<b>%{text}</b><br><br>" + "\u0394\u0394G_Exp: %{x:.2f}<br>" + "\u0394\u0394G_Pred: %{y:.2f}<br>"
            name_plot = '\u0394\u0394G'
        else:
            hov_label = "<b>%{text}</b><br><br>" + "\u0394G_Exp: %{x:.2f}<br>" + "\u0394G_Pred: %{y:.2f}<br>"
            name_plot = '\u0394G'

        figure.add_trace(go.Scatter(x=x, y=y,
                                    text=hover_text,
                                    hovertemplate=hov_label,
                                    name=name_plot,
                                    mode='markers',
                                    xaxis='x1',
                                    yaxis='y1',
                                    showlegend=False,
                                    error_x=dict(
                                        type='data',
                                        array=err_x,
                                        color='blue',
                                        visible=True),
                                    error_y=dict(
                                        type='data',
                                        array=err_y,
                                        color='blue',
                                        visible=True),
                                    marker=dict(color='blue', size=8)
                                    ),
                         row=1, col=col)

        return

    #######################

    skip_plot_methods = ['RMSE', 'RRMSE', 'RHO']

    display_units = units

    raw_results = []

    if len(exp_data_dic) != len(predicted_data_dic):
        raise ValueError("Experimental vs Predicted data mismatch number: {} vs {}".
                         format(len(exp_data_dic), len(predicted_data_dic)))

    for edge_name, exp_val in exp_data_dic.items():
        ligA_name = edge_name.split()[0]
        ligB_name = edge_name.split()[2]

        res = wrangle.Result(ligA_name, ligB_name,
                             exp_val[0], exp_val[1],
                             predicted_data_dic[edge_name][0], predicted_data_dic[edge_name][1], 0.0)

        raw_results.append(res)

    network = wrangle.FEMap(raw_results)

    network.check_weakly_connected()

    # Relative Binding Affinity
    exp_DDG = np.asanyarray([x.exp_DDG for x in network.results])
    dexp_DDG = np.asanyarray([x.dexp_DDG for x in network.results])
    ligA_names = [x.ligandA for x in network.results]

    pred_DDG = np.asanyarray([y.calc_DDG for y in network.results])
    dpred_DDG = np.asanyarray([y.mbar_dDDG for y in network.results])
    ligB_names = [y.ligandB for y in network.results]

    hov_edge = []
    for lna, lnb in zip(ligA_names, ligB_names):
        hov_edge.append("{} to {}".format(lna, lnb))

    if DDG_symmetrize:
        exp_DDG = np.asarray([x.exp_DDG for x in network.results] + [-x.exp_DDG for x in network.results])
        dexp_DDG = np.asarray([x.dexp_DDG for x in network.results] + [x.dexp_DDG for x in network.results])

        pred_DDG = np.asarray([y.calc_DDG for y in network.results] + [-y.calc_DDG for y in network.results])
        dpred_DDG = np.asarray([y.mbar_dDDG for y in network.results] + [y.mbar_dDDG for y in network.results])
        for lna, lnb in zip(ligB_names, ligA_names):
            hov_edge.append("{} to {}".format(lna, lnb))

    # Relative statistics
    relative_statistics = calculate_statistics(exp_DDG, pred_DDG, plot_type="ddG")

    if network.weakly_connected:
        figure = make_subplots(rows=2, cols=2,
                               subplot_titles=("Relative Binding Affinity",
                                               "Centered Binding Affinity"),
                               vertical_spacing=0.0,
                               row_heights=[0.5, 0.5],
                               )
    else:
        figure = make_subplots(rows=2, cols=2, subplot_titles=("Relative Binding Affinity",
                                                               "Centered Binding Affinity Not Available"))

    # Axis range:
    conc_DDG = np.concatenate((exp_DDG, pred_DDG))

    min_range_DDG = np.min(conc_DDG)
    max_range_DDG = np.max(conc_DDG)

    conc_dDDG = np.concatenate((dexp_DDG, dpred_DDG))

    max_dDDG = np.max(conc_dDDG)

    # Correct for the error bars
    min_range_DDG -= max_dDDG
    max_range_DDG += max_dDDG

    padding = int(round((np.abs(min_range_DDG) + np.abs(max_range_DDG)) * 0.05))

    min_range_DDG -= padding
    max_range_DDG += padding

    range_axis_DDG = [min_range_DDG, max_range_DDG]

    sub_plot(exp_DDG, dexp_DDG, pred_DDG, dpred_DDG, figure, hov_edge, 1, range_axis_DDG)

    xaxis1 = dict(
        zeroline=True,
        showgrid=True,
        tickfont=dict(
            size=18,
        ),
        range=range_axis_DDG,
        constrain='domain',
        title="Experimental \u0394\u0394G kJ/mol",
        titlefont=dict(
            size=20
        )
    )

    yaxis1 = dict(
        zeroline=True,
        showgrid=True,
        tickfont=dict(
            size=18,
        ),
        scaleanchor="x1",
        scaleratio=1,
        range=range_axis_DDG,
        title="Predicted \u0394\u0394G kJ/mol",
        titlefont=dict(
            size=20
        )
    )

    figure.update_layout(
        title='Method used to estimate \u0394\u0394G: {}'.format(method),
        xaxis1=xaxis1,
        yaxis1=yaxis1,
        autosize=False,
        width=1100,
        height=1000,
    )

    if network.weakly_connected:

        affinity_dic = dict()

        # Absolute Binding Affinity from Hanna's code
        graph = network.graph

        exp_DG = np.asarray([node[1]['f_i_exp'] for node in graph.nodes(data=True)])
        pred_DG = np.asarray([node[1]['f_i_calc'] for node in graph.nodes(data=True)])
        dexp_DG = np.asarray([node[1]['df_i_exp'] for node in graph.nodes(data=True)])
        dpred_DG = np.asarray([node[1]['df_i_calc'] for node in graph.nodes(data=True)])

        # centralising
        # this should be replaced by providing one experimental result
        exp_DG = exp_DG - np.mean(exp_DG)
        pred_DG = pred_DG - np.mean(pred_DG)

        if relative_statistics is not None:
            # Absolute statistics
            absolute_statistics = calculate_statistics(exp_DG, pred_DG, plot_type="dG")

        # Assuming that dic is in Order
        hov_ligs = []
        for name, id in network._name_to_id.items():
            hov_ligs.append(name)

        if not (len(hov_ligs) == len(exp_DG) == len(dexp_DG) == len(pred_DG) == len(dpred_DG)):
            raise ValueError("List size mismatch")

        for name, exp, exp_err, pred, pred_err in zip(hov_ligs, exp_DG, dexp_DG, pred_DG, dpred_DG):
            affinity_dic[name] = [exp, exp_err, pred, pred_err]

        # for name, l in affinity_dic.items():
        #     print("name {} exp {} +- {} pred {} +- {}".format(name, l[0], l[1], l[2], l[3]))

        # Axis range:
        conc_DG = np.concatenate((exp_DG, pred_DG))

        min_range_DG = np.min(conc_DG)
        max_range_DG = np.max(conc_DG)

        conc_dDG = np.concatenate((dexp_DG, dpred_DG))
        max_dDG = np.max(conc_dDG)

        # Correct for the error bars
        min_range_DG -= max_dDG
        max_range_DG += max_dDG

        padding = int(round((np.abs(min_range_DG) + np.abs(max_range_DG)) * 0.05))

        min_range_DG -= padding
        max_range_DG += padding

        range_axis_DG = [min_range_DG, max_range_DG]

        sub_plot(exp_DG, dexp_DG, pred_DG, dpred_DG, figure, hov_ligs, 2, range_axis_DG)

        xaxis2 = dict(
            zeroline=True,
            showgrid=True,
            tickfont=dict(
                size=18,
            ),

            range=range_axis_DG,
            constrain='domain',
            title="Experimental \u0394G kJ/mol",
            titlefont=dict(
                size=20
            )
        )

        yaxis2 = dict(
            zeroline=True,
            showgrid=True,
            tickfont=dict(
                size=18,
            ),
            scaleanchor="x2",
            scaleratio=1,
            range=range_axis_DG,
            title="Predicted \u0394G kJ/mol",
            titlefont=dict(
                size=20
            )
        )

        figure.update_layout(
            xaxis2=xaxis2,
            yaxis2=yaxis2
        )
    else:
        print('Graph is not connected enough to compute absolute values')

    data = []
    methods = []
    values = []

    if relative_statistics is not None:
        for meth, val in relative_statistics.items():

            # Skip these methods
            if meth in skip_plot_methods:
                continue
            if meth == 'RAE':
                meth = "Relative MAE"
            if meth == 'RRMSE':
                meth = "Relative RMSE"
            if meth == 'R2':
                meth = "Pearson's R\N{SUPERSCRIPT TWO}"
            if meth == 'RHO':
                meth = "Spearman's \u03C1"
            if meth == "KTAU":
                meth = "Kendall's \u03C4"

            methods.append(meth)

            vl_line = "{:.2f} [ {:.2f} {:.2f}]".format(val['data'], val['low'], val['high'])
            values.append(vl_line)

        data.append(methods)
        data.append(values)

        data.append([])

        values = []

        if network.weakly_connected:
            for meth, val in absolute_statistics.items():

                if meth in skip_plot_methods:
                    continue

                vl_line = "{:.2f} [ {:.2f} {:.2f}]".format(val['data'], val['low'], val['high'])
                values.append(vl_line)

        data.append(methods)
        data.append(values)

        # Add Table with statistics. Data in kJ/mol
        alpha_color = 'rgba(127, 166, 238, 0.4)'
        figure.add_trace(
            go.Table(
                header=dict(
                    values=["Method", "\u0394\u0394G kJ/mol", "", "Method", "\u0394G kJ/mol"],
                    font=dict(size=14),
                    align="center",
                    fill_color=alpha_color,
                ),
                cells=dict(
                    values=data,
                    align="center"),

                domain=dict(x=[0, 1],
                            y=[0, 0.35])
            )
        )

    # figure.show()

    with TemporaryDirectory() as outdir:

        fn = os.path.join(outdir, "report.html")

        plotly.offline.plot(figure, filename=fn)

        with open(fn, 'r') as f:
            report_string_list = f.readlines()

        # Add the following lines for Orion CSS Plotly integration
        report_str = ""

        for idx in range(0, len(report_string_list)):
            line = report_string_list[idx]
            if "<head>" in line:
                report_str += """<head><meta charset="utf-8" />   
                   <style>
                       .modebar { display:unset; }
                       .js-plotly-plot { min-height: 900px; }
                       .js-plotly-plot { min-width: 900px; }
                   </style></head>\n"""

            report_str += line

    if network.weakly_connected:
        return report_str, affinity_dic
    else:
        return report_str, None


def upload_gmx_files(tar_fn, record, shard_name=""):

    if not os.path.isfile(tar_fn):
        raise ValueError("The filename does not exist: {}".format(tar_fn))

    if in_orion():

        if not record.has_field(Fields.collections):
            raise ValueError("Missing Collection Field")

        # session = APISession
        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        collections_dic = record.get_value(Fields.collections)

        nes_collection_id = collections_dic[CollectionsNames.nes]

        collection = session.get_resource(ShardCollection, nes_collection_id)

        shard = try_hard_to_create_shard(collection, tar_fn, name=shard_name)

        record.set_value(Fields.extra_data_tar, shard.id)

        shard.close()

    else:
        record.set_value(Fields.extra_data_tar, tar_fn)

    return True


def download_gmx_file(mdrecord):

    if not mdrecord.has_field(Fields.extra_data_tar):
        raise ValueError("Extra data tar file has not been found on the record")

    extra_data = mdrecord.get_value(Fields.extra_data_tar)

    if in_orion():

        if not mdrecord.has_field(Fields.collections):
            raise ValueError("Collection field has not been found on the record")

        collection_dic = mdrecord.get_value(Fields.collections)
        nes_collection_id = collection_dic[CollectionsNames.nes]

        # session = APISession

        session = OrionSession(
            requests_session=get_session(
                retry_dict={
                    403: 5,
                    404: 20,
                    409: 45,
                    460: 15,
                    500: 2,
                    502: 45,
                    503: 45,
                    504: 45,
                }
            )
        )

        collection = session.get_resource(ShardCollection, nes_collection_id)

        shard = session.get_resource(Shard(collection=collection), extra_data)

        fn = os.path.join(mdrecord.cwd, 'extra_data.tar')

        try_hard_to_download_shard(shard, fn)

        shard.close()

        extra_data = fn

    return extra_data
