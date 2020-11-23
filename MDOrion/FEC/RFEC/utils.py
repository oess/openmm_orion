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

from simtk import openmm

import glob

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


def extract_trj_velocities(trj_fn,  trj_type="OpenMM"):

    if trj_type == "OpenMM":

        f = h5py.File(trj_fn, 'r')

        # velocities = f.get('velocities').value
        velocities = f['velocities'][()]

        f.close()

    # TODO MISSING GMX PART

    return velocities


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

    stage = mdrecord.get_stage_by_name('last')
    trj_field = stage.get_field(Fields.trajectory.get_name())
    trj_meta = trj_field.get_meta()
    md_engine = trj_meta.get_attribute(Meta.Annotation.Description)

    # Velocities array shape = (frame, atoms, 3) units nm/ps
    # TODO MISSING GMX PART FOR EQUILIBRIUM SIMULATION
    velocities = extract_trj_velocities(trj_fn, trj_type=md_engine)

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
                                         gen_vel='no',
                                         continue_sim='no',
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


def plot_work_pdf(f_bound, r_bound, f_unbound, r_unbound, results, title, edge_depiction_image):

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
        x_left = 'x'+str(row+1)
        x_right = 'x'+str(row+2)
        y_left = 'y'+str(row+1)
        y_right = 'y'+str(row+2)

        legend = False

        if row == 2:
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

    fig = make_subplots(rows=4, cols=2,
                        horizontal_spacing=0.01,
                        vertical_spacing=0.1,
                        column_widths=[0.7, 0.3],
                        shared_yaxes=True)
    # Bound Plot
    make_plots(bound_frames_f, bound_frames_r, bound_work_f, bound_work_r, fig, row=2)

    # Unbound Plot
    make_plots(unbound_frames_f, unbound_frames_r, unbound_work_f, unbound_work_r, fig, row=3)

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
    #
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
                        y=[0, 0.15])
        )

    )

    # Add edge depiction
    with TemporaryDirectory() as outdir:

        fn = os.path.join(outdir, "depiction.png")

        oedepict.OEWriteImage(fn, edge_depiction_image)

        image1 = Image.open(fn)

        fig.add_layout_image(
            dict(
                source=image1,
                xref="x1", yref="y1",
                x=-0.5, y=5.8,
                sizex=10, sizey=8,
                xanchor="left", yanchor="top",
                layer="above", opacity=1), row=1, col=1

        )

    # Muatuation Axes
    xaxis1 = dict(
        zeroline=False,
        showgrid=False,
        title="Mutation",
        showticklabels=False,
        visible=True

    )

    yaxis1 = dict(
        zeroline=False,
        showgrid=False,
        title="Molecule",
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
        title="Work Bound kJ/mol"
    )

    # Bound data PDF
    xaxis4 = dict(
        zeroline=False,
        showgrid=False,
        visible=True,
        title="PDF"
    )

    yaxis4 = dict(
        zeroline=False,
        showgrid=False,
        visible=False
    )

    # Bound data Frames vs Work
    xaxis5 = dict(
        zeroline=True,
        showgrid=True,
        title="Frames"
    )

    yaxis5 = dict(
        zeroline=True,
        showgrid=True,
        title="Work Unbound kJ/mol"
    )

    # Bound data PDF
    xaxis6 = dict(
        zeroline=False,
        showgrid=False,
        title="PDF"
    )

    yaxis6 = dict(
        zeroline=False,
        showgrid=False,
        visible=False
    )

    fig.update_layout(

        legend=dict(x=0.7, y=0.8),

        barmode='overlay',
        hovermode="closest",
        bargap=0,
        title="Non Equilibrium Switching Results: {}".format(title),

        # Mutation axes
        xaxis1=xaxis1,
        yaxis1=yaxis1,

        # Bound data frames
        xaxis3=xaxis3,
        yaxis3=yaxis3,

        # Bound data PDF
        xaxis4=xaxis4,
        yaxis4=yaxis4,

        # Unbound data frames
        xaxis5=xaxis5,
        yaxis5=yaxis5,

        # Unbound data PDF
        xaxis6=xaxis6,
        yaxis6=yaxis6,

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

    return report_str
