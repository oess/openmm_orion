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

from simtk import unit

from simtk.openmm import app

import simtk

from MDOrion.MDEngines.utils import (MDSimulations,
                                     md_keys_converter)

import numpy as np

from MDOrion.MDEngines.Gromacs.gromacs_templates import (gromacs_minimization,
                                                         gromacs_nvt_npt,
                                                         gromacs_pos_restraints)
import subprocess

import parmed

import copy

from oeommtools import utils as oeommutils

import tarfile

import os

import io

from itertools import accumulate

import itertools

from collections import OrderedDict

from MDOrion.Standards import MDEngines

from subprocess import STDOUT, PIPE, Popen, DEVNULL

import math

from parmed.gromacs.gromacstop import _Defaults


class GromacsSimulations(MDSimulations):

    def __init__(self, mdstate, parmed_structure, opt):
        super().__init__(mdstate, parmed_structure, opt)

        velocities = mdstate.get_velocities()
        box = mdstate.get_box_vectors()

        if box is not None:
            omm_system = parmed_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                                       nonbondedCutoff=10.0 * unit.angstroms,
                                                       constraints=None,
                                                       removeCMMotion=False,
                                                       rigidWater=False)
        else:
            omm_system = parmed_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                       constraints=None,
                                                       removeCMMotion=False,
                                                       rigidWater=False)
        # Define unique atom types
        atom_types_dic = {}
        count_id = 0

        # # Copy the topology and positions
        topology = parmed_structure.topology
        positions = parmed_structure.positions

        # Parmed GMX set 1-4 scaling
        defaults = _Defaults(fudgeLJ=0.5, fudgeQQ=0.8333, gen_pairs='yes')
        parmed_structure.defaults = defaults

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

        # Rename all the water molecules to avoid Gromacs
        # settling errors
        for r in topology.residues():
            if check_water(r):
                h1 = False
                for a in r.atoms():
                    if a.element.atomic_number == 1:
                        if not h1:
                            a.name = 'H1'
                            h1 = True
                        else:
                            a.name = 'H2'
                    else:
                        a.name = 'O'

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
                        else:
                            a.id = 'O' + str(count_id)
                            count_id += 1
                            atom_types_dic[r.name + a.name] = a.id

        # Define a new parmed structure with the new unique atom types
        new_system_structure = parmed.openmm.load_topology(topology,
                                                           system=omm_system,
                                                           xyz=positions)

        # Re-set positions, velocities and box
        new_system_structure.positions = parmed_structure.positions
        new_system_structure.velocities = parmed_structure.velocities
        new_system_structure.box = parmed_structure.box
        new_system_structure.defaults = defaults

        self.stepLen = 0.002 * unit.picoseconds

        opt['timestep'] = self.stepLen

        if box is not None:

            box_v = parmed_structure.box_vectors.value_in_unit(unit.angstrom)
            box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])

            min_box = np.min(box_v)

            threshold = (min_box / 2.0) * 0.85

            if opt['nonbondedCutoff'] < threshold:
                cutoff_distance = opt['nonbondedCutoff'] * unit.angstroms
            else:
                opt['Logger'].warn("[{}] Cutoff Distance too large for the box size. Set the cutoff distance "
                                   "to {} A".format(opt['CubeTitle'], threshold))

                cutoff_distance = threshold * unit.angstroms

            cutoff = cutoff_distance

        # Centering the system
        if opt['center'] and box is not None:
            opt['Logger'].info("[{}] Centering is On".format(opt['CubeTitle']))
            # Numpy array in A
            coords = new_system_structure.coordinates
            # System Center of Geometry
            cog = np.mean(coords, axis=0)
            # System box vectors
            box_v = new_system_structure.box_vectors.in_units_of(unit.angstrom) / unit.angstrom
            box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])
            # Translation vector
            delta = box_v / 2 - cog
            # New Coordinates
            new_coords = coords + delta
            new_system_structure.coordinates = new_coords
            positions = new_system_structure.positions
            mdstate.set_positions(positions)

        if box is not None:
            pbc = 'xyz'
            nslist = 10
        else:
            pbc = 'no'
            cutoff = 0.0 * unit.angstrom
            nslist = 0

        if opt['SimType'] == 'min':

            if opt['steps'] == 0:
                max_minimization_steps = 100000
            else:
                max_minimization_steps = opt['steps']

            mdp_template = gromacs_minimization.format(
                nsteps=max_minimization_steps,
                nslist=1,
                cutoff=cutoff.value_in_unit(unit.nanometer),
                pbc=pbc,
            )

        if opt['SimType'] in ['nvt', 'npt']:

            if velocities is None:
                opt['Logger'].info('[{}] GENERATING a new starting State'.format(opt['CubeTitle']))
                gen_vel = 'yes'
            else:
                gen_vel = 'no'

            if opt['SimType'] == "npt":
                # If restraints
                if opt['restraints']:
                    pcoupl = 'berendsen'
                else:
                    pcoupl = 'Parrinello-Rahman'
            else:  # nvt ensemble do not use any pressure
                pcoupl = 'no'
                # This is not used
                opt['pressure'] = 0.0

            if opt['reporter_interval']:
                reporter_steps = int(round(opt['reporter_interval'] / (
                        opt['timestep'].in_units_of(unit.nanoseconds) / unit.nanoseconds)))
            else:
                reporter_steps = 0

            # Convert simulation time in steps
            opt['steps'] = int(round(opt['time'] / (self.stepLen.in_units_of(unit.nanoseconds) / unit.nanoseconds)))

            if opt['trajectory_interval']:
                trajectory_steps = int(round(opt['trajectory_interval'] / (
                        opt['timestep'].in_units_of(unit.nanoseconds) / unit.nanoseconds)))

            elif opt['trajectory_frames']:

                if opt['steps'] < opt['trajectory_frames']:
                    raise ValueError(" The selected number of frames {} will exceed the total produced md steps {}".
                                     format(opt['trajectory_frames'], opt['steps']))

                trajectory_steps = int(math.floor(opt['steps'] / opt['trajectory_frames']))
            else:
                trajectory_steps = 0

            # Constraints
            constraints = md_keys_converter[MDEngines.Gromacs]['constraints'][opt['constraints']]

            mdp_template = gromacs_nvt_npt.format(
                nsteps=opt['steps'],
                timestep=self.stepLen.in_units_of(unit.picoseconds) / unit.picoseconds,
                reporter_steps=reporter_steps,
                trajectory_steps=trajectory_steps,
                constraints=constraints,
                nslist=nslist,
                cutoff=cutoff.in_units_of(unit.nanometer) / unit.nanometer,
                temperature=opt['temperature'],
                pcoupl=pcoupl,
                pressure=opt['pressure'],
                gen_vel=gen_vel,
                pbc=pbc
            )

        opt['Logger'].info("Output Directory {}".format(opt['out_directory']))

        opt['outfname'] = opt['system_title'] + '_' + str(opt['system_id']) + '-' + opt['suffix']

        # Gromacs file names
        opt['grm_top_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".top")

        opt['grm_top_res_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + "_res.top")
        opt['grm_ndx_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".ndx")

        opt['grm_gro_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".gro")
        opt['grm_pdb_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".pdb")
        opt['grm_tpr_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".tpr")
        opt['mdp_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".mdp")
        opt['grm_def_fn'] = os.path.join(opt['out_directory'], opt['outfname']+"_run")
        opt['grm_log_fn'] = opt['grm_def_fn'] + '.log'
        opt['grm_trj_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".trr")
        # opt['grm_trj_comp_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".xtc")

        opt['mdp_template'] = mdp_template

        # Generate coordinate file
        new_system_structure.save(opt['grm_gro_fn'], overwrite=True)
        new_system_structure.save(opt['grm_pdb_fn'], overwrite=True)

        # Generate Gromacs .mdp configuration files
        with open(opt['mdp_fn'], 'w') as of:
            of.write(mdp_template)

        apply_restraints = False

        # Select atom to restraint
        if opt['restraints']:

            res_atom_list = sorted(oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['restraints']))

            if res_atom_list:

                apply_restraints = True

                res_atom_set = set(res_atom_list)

                pmd_split = new_system_structure.split()

                MAX_DIGITS = len(str(len(new_system_structure.atoms)))

                atom_count_list = []

                for tp in pmd_split:
                    for str_idx in tp[1]:
                        atom_count_list.insert(str_idx, len(tp[0].atoms))

                cumul_list = list(accumulate(atom_count_list))

                systems = []

                for tp in pmd_split:
                    tmp_dic = OrderedDict()
                    for str_idx in sorted(tp[1]):
                        if str_idx == 0:
                            tmp_dic[0] = {x for x in range(0, atom_count_list[0])}
                        else:
                            tmp_dic[str_idx] = {x for x in range(cumul_list[str_idx - 1], cumul_list[str_idx])}

                    systems.append(tmp_dic)

        # Apply restraints
        if apply_restraints:

            if opt['restraint_to_reference']:
                opt['Logger'].info("[{}] Restraint to the Reference State Enabled".format(opt['CubeTitle']))
                reference_positions = opt['reference_state'].get_positions()
                coords = np.array(reference_positions.value_in_unit(unit.angstrom))
                # System Center of Geometry
                cog = np.mean(coords, axis=0)

                # System box vectors
                box_v = opt['reference_state'].get_box_vectors().value_in_unit(unit.angstrom)
                box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])

                # Translation vector
                delta = box_v / 2 - cog
                # New Coordinates
                corrected_reference_positions = coords + delta

                tmp_pmd_structure = new_system_structure.__copy__()

                tmp_pmd_structure.coordinates = corrected_reference_positions
                opt['grm_gro_ref_restraint_fn'] = os.path.join(opt['out_directory'],
                                                               opt['outfname'] + "_ref_restraint.gro")
                tmp_pmd_structure.save(opt['grm_gro_ref_restraint_fn'], overwrite=True)

            # Restraints weight
            res_wgt = opt['restraintWt'] * unit.kilocalories_per_mole / (unit.angstroms ** 2)

            # grm unit
            res_wgt_grm = res_wgt.value_in_unit(unit.kilojoule_per_mole / (unit.nanometer ** 2))

            system_order = dict()

            restraint_dic = dict()

            for idx in range(0, len(systems)):

                # print("IDX = {}".format(idx))

                s_union = set.union(*[s for s in systems[idx].values()])

                # for k, v in systems[idx].items():
                #     pmd_structures[k] = new_fix_atom_type_pmd[list(v)]

                # The intersection between the atom_index set and
                # the restrained atom index set is empty
                if not s_union.intersection(res_atom_set):
                    # No restraint file needs to be generated in this case

                    for k, v in systems[idx].items():
                        system_order[k] = 'PMD' + '_' + str(idx)
                    continue

                # The intersection between the atom_index set and
                # the restrained atom index set is equal to the atom_index set
                # In this case all the atom_index set needs to restrained
                # generating one restraint file for the atom_index set part
                elif s_union.intersection(res_atom_set) == s_union:

                    sn = 'PMD' + '_' + str(idx)

                    for k, v in systems[idx].items():
                        system_order[k] = sn

                    sfn = os.path.join(opt['out_directory'], "system" + str(idx) + '.itp')

                    with open(sfn, 'w') as of:
                        of.write("; FULL\n")
                        of.write("[ position_restraints ]\n")
                        of.write("; ai  funct  fcx    fcy    fcz\n")
                        of.write(gromacs_pos_restraints.format(1, 1,
                                                               res_wgt_grm,
                                                               res_wgt_grm,
                                                               res_wgt_grm,
                                                               digits=MAX_DIGITS))

                        restraint_dic[sn] = sfn

                else:

                    for k, v in systems[idx].items():

                        if not v.intersection(res_atom_set):
                            # print("EMPTY")
                            sn = 'PMD' + '_' + str(idx)
                            system_order[k] = sn
                            continue

                        # The intersection is total. Generates a restraint file with the
                        # just one line that restraints all the part
                        elif v.intersection(res_atom_set) == v:
                            # print("FULL")
                            sn = 'PMD' + '_' + str(idx) + '_' + str(k)

                            system_order[k] = sn

                            sfn = os.path.join(opt['out_directory'], "system" + str(idx) + '_' + str(k) + '.itp')

                            with open(sfn, 'w') as of:
                                of.write("; FULL\n")
                                of.write("[ position_restraints ]\n")
                                of.write("; ai  funct  fcx    fcy    fcz\n")
                                of.write(gromacs_pos_restraints.format(1, 1,
                                                                       res_wgt_grm,
                                                                       res_wgt_grm,
                                                                       res_wgt_grm,
                                                                       digits=MAX_DIGITS))

                            restraint_dic[sn] = sfn

                        # The intersection is partial. Generates a restraint file with the
                        # atom index to restraint
                        else:  # Atom Indexes are fixed after
                            # print("PARTIAL")

                            sn = 'PMD' + '_' + str(idx) + '_' + str(k)

                            system_order[k] = sn

                            # Intersection
                            intersec = v.intersection(res_atom_set)

                            sfn = os.path.join(opt['out_directory'], "system" + str(idx) + '_' + str(k) + '.itp')

                            with open(sfn, 'w') as of:
                                of.write("; PARTIAL\n")
                                of.write("[ position_restraints ]\n")
                                of.write("; ai  funct  fcx    fcy    fcz\n")
                                for ai in list(sorted(intersec)):
                                    of.write(gromacs_pos_restraints.format(ai + 1, 1,
                                                                           res_wgt_grm,
                                                                           res_wgt_grm,
                                                                           res_wgt_grm,
                                                                           digits=MAX_DIGITS))

                            restraint_dic[sn] = sfn

            list_all = [system_order[k] for k in sorted(system_order)]
            # print(list_all)
            occ = [(k, len((list(g)))) for k, g in itertools.groupby(list_all)]
            # print(occ)

            # Generate topology files
            new_system_structure.save(opt['grm_top_fn'], overwrite=True)

            f = open(opt['grm_top_fn'], 'r')

            header_str = ""
            for line in f:
                if '[ moleculetype ]' in line:
                    break
                header_str += line

            f.close()

            header_list = header_str.split("\n")

            for idx in range(0, len(header_list)):
                if '[ defaults ]' in header_list[idx]:
                    header_list[idx + 2] = '1               2               yes             0.5          0.83333333  \n'
                    break

            header_str = "\n".join(header_list)

            unique_pmd_names = list(dict.fromkeys([pair[0] for pair in occ]))

            # print(unique_pmd_names)

            restraint_str = """

            ; Include Position restraint file
            #ifdef POSRES
            #include "{}"
            #endif

            """

            moltype_str = ""

            for pmd_name in unique_pmd_names:
                pmd_id = int(pmd_name.split('_')[1])

                fn = os.path.join(opt['out_directory'], pmd_name + '.top')

                pmd_split[pmd_id][0].save(fn, overwrite=True)

                f = open(fn, 'r')

                tmp_str_list = f.readlines()
                f.close()

                capture = False
                block = []
                for ln in tmp_str_list:
                    if '[ moleculetype ]' in ln:
                        capture = True
                    if '[ system ]' in ln:
                        capture = False
                    if capture:
                        block.append(ln)

                tmp_str_exc = block[2].split()[-1]
                block[2] = pmd_name + '             {}\n'.format(tmp_str_exc)
                moltype_str += ''.join(block)

                if pmd_name in restraint_dic:
                    moltype_str += restraint_str.format(restraint_dic[pmd_name])

            system_header_str = """
            [ system ]
            ; Name
            Generic title

            [ molecules ]
            ; Compound       #mols

            """

            system_count_str = ""
            curr_total_atoms = 0

            for pair in occ:
                system_count_str += pair[0] + "          {}\n".format(pair[1])

                # Fix Atom index in the restraint files
                if pair[0] in restraint_dic:

                    fn = os.path.join('opt_directory', restraint_dic[pair[0]])

                    f = open(fn, 'r')

                    restraint_str = ""

                    for ln in f:
                        if "FULL" in ln:
                            break
                        elif ";" in ln:
                            restraint_str += ln
                        elif "[" in ln:
                            restraint_str += ln
                        else:
                            list_ln = ln.split()
                            list_ln[0] = str(int(list_ln[0]) - curr_total_atoms)

                            restraint_str += gromacs_pos_restraints.format(int(list_ln[0]),
                                                                           int(list_ln[1]),
                                                                           float(list_ln[2]),
                                                                           float(list_ln[3]),
                                                                           float(list_ln[4]),
                                                                           digits=MAX_DIGITS)
                    f.close()
                    if restraint_str:
                        fn = os.path.join('opt_directory', restraint_dic[pair[0]])
                        f = open(fn, 'w')
                        f.write(restraint_str)
                        f.close()

                pmd_id = int(pair[0].split('_')[1])

                curr_total_atoms += int(pair[1]) * len(pmd_split[pmd_id][0].atoms)

            system_str = system_header_str + system_count_str

            f = open(opt['grm_top_res_fn'], 'w')

            f.write(header_str)
            f.write(moltype_str)
            f.write(system_str)
            f.close()

            # Generate Atom index file .ndx used to apply the restraints
            chunk = 15
            count = 0
            f = open(opt['grm_ndx_fn'], "a+")

            atom_idx = range(1, len(new_system_structure.atoms) + 1)
            f.write("[ System ]\n")

            for idx in atom_idx:
                f.write("{:>{digits}}".format(idx, digits=MAX_DIGITS))
                if count == chunk - 1:
                    count = 0
                    f.write("\n")
                else:
                    f.write(" ")
                    count += 1
            f.write("\n")
            f.close()

            for k, v in restraint_dic.items():

                fn = os.path.join(opt['out_directory'], v)
                f = open(fn, 'r')

                atom_idx = []
                for l in f:
                    if "FULL" in l:
                        pmd_id = int(k.split('_')[1])
                        atom_idx = range(1, len(pmd_split[pmd_id][0].atoms) + 1)
                        break
                    elif ";" in l:
                        continue
                    elif "[" in l:
                        continue
                    else:
                        list_ln = l.split()
                        atom_idx.append(list_ln[0])
                f.close()
                f = open(opt['grm_ndx_fn'], "a+")
                f.write("[ {} ]\n".format(k))

                count = 0
                for idx in atom_idx:
                    f.write("{:>{digits}}".format(idx, digits=MAX_DIGITS))
                    if count == chunk - 1:
                        count = 0
                        f.write("\n")
                    else:
                        f.write(" ")
                        count += 1
                f.write("\n")
                f.close()

            opt['Logger'].info("[{}] RESTRAINT mask applied to: {}"
                               "\tRestraint weight: {}".format(opt['CubeTitle'],
                                                               opt['restraints'],
                                                               opt[
                                                                   'restraintWt'] * unit.kilocalories_per_mole / unit.angstroms ** 2))

            # Select atom to restraint
            opt['Logger'].info("[{}] Number of restraint atoms: {}".format(opt['CubeTitle'], len(res_atom_list)))

        # Generate Gromacs .tpr file
        if apply_restraints:

            if opt['restraint_to_reference']:
                subprocess.check_call(['gmx',
                                       'grompp',
                                       '-f', opt['mdp_fn'],
                                       '-c', opt['grm_gro_fn'],
                                       '-r', opt['grm_gro_ref_restraint_fn'],
                                       '-p', opt['grm_top_res_fn'],
                                       '-n', opt['grm_ndx_fn'],
                                       '-o', opt['grm_tpr_fn'],
                                       '-maxwarn', b'5'])
            else:
                subprocess.check_call(['gmx',
                                       'grompp',
                                       '-f', opt['mdp_fn'],
                                       '-c', opt['grm_gro_fn'],
                                       '-r', opt['grm_gro_fn'],
                                       '-p', opt['grm_top_res_fn'],
                                       '-n', opt['grm_ndx_fn'],
                                       '-o', opt['grm_tpr_fn'],
                                       '-maxwarn', b'5'])

        else:
            # Generate topology files
            new_system_structure.save(opt['grm_top_fn'], overwrite=True)

            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', opt['mdp_fn'],
                                   '-c', opt['grm_gro_fn'],
                                   '-p', opt['grm_top_fn'],
                                   '-o', opt['grm_tpr_fn'],
                                   '-maxwarn', b'5'])

        self.mdstate = mdstate
        self.opt = opt

        return

    def run(self):

        # Run Gromacs
        if self.opt['verbose']:

            subprocess.check_call(['gmx',
                                   'mdrun',
                                   '-v',
                                   '-s', self.opt['grm_tpr_fn'],
                                   '-deffnm', self.opt['grm_def_fn'],
                                   '-o', self.opt['grm_trj_fn']
                                   ])
        else:
            p = Popen(['gmx',
                       'mdrun',
                       '-v',
                       '-s', self.opt['grm_tpr_fn'],
                       '-deffnm', self.opt['grm_def_fn'],
                       '-o', self.opt['grm_trj_fn']
                       ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)

            p.communicate()

        if self.opt['SimType'] in ['nvt', 'npt']:

            if self.opt['reporter_interval']:

                with(io.open(self.opt['grm_log_fn'], 'r', encoding='utf8', errors='ignore')) as fr:
                    log_string = fr.read()

                self.opt['str_logger'] += '\n'+log_string

            # Save trajectory files
            if self.opt['trajectory_interval'] or self.opt['trajectory_frames']:

                # Generate whole system trajectory
                # p = subprocess.Popen(['gmx',
                #                       'trjconv',
                #                       '-f', self.opt['grm_trj_fn'],
                #                       '-s', self.opt['grm_tpr_fn'],
                #                       '-o', self.opt['grm_trj_comp_fn'],
                #                       '-pbc', b'whole'],
                #                      stdin=subprocess.PIPE)
                #
                # # Select the entire System
                # p.communicate(b'0')

                # Tar the files dir with its content:
                tar_fn = self.opt['trj_fn']

                with tarfile.open(tar_fn, mode='w:gz') as archive:

                    archive.add(self.opt['grm_gro_fn'], arcname=os.path.basename(self.opt['grm_gro_fn']))
                    archive.add(self.opt['grm_pdb_fn'], arcname=os.path.basename(self.opt['grm_pdb_fn']))
                    archive.add(self.opt['grm_top_fn'], arcname=os.path.basename(self.opt['grm_top_fn']))
                    archive.add(self.opt['grm_trj_fn'], arcname=os.path.basename(self.opt['grm_trj_fn']))
                    # archive.add(self.opt['grm_trj_comp_fn'], arcname=os.path.basename(self.opt['grm_trj_comp_fn']))
                    archive.add(self.opt['grm_log_fn'], arcname=os.path.basename(self.opt['grm_log_fn']))
                    archive.add(self.opt['grm_tpr_fn'], arcname=os.path.basename(self.opt['grm_tpr_fn']))
                    archive.add(self.opt['mdp_fn'], arcname=os.path.basename(self.opt['mdp_fn']))

                # with tarfile.open(tar_fn, mode='w:gz') as archive:
                #     archive.add(self.opt['out_directory'], arcname=os.path.basename(self.opt['out_directory']))

        return

    def update_state(self):

        # top_fn = self.opt['grm_top_fn']

        gro_fn = os.path.join(self.opt['out_directory'], self.opt['grm_def_fn'] + '.gro')

        gro_structure = parmed.gromacs.GromacsGroFile().parse(gro_fn, skip_bonds=True)

        new_mdstate = copy.deepcopy(self.mdstate)

        new_mdstate.set_positions(gro_structure.positions)

        if gro_structure.box_vectors is not None:
            new_mdstate.set_box_vectors(gro_structure.box_vectors)

        if self.opt['SimType'] in ['nvt', 'npt']:
            new_mdstate.set_velocities(gro_structure.velocities * simtk.unit.angstrom/simtk.unit.picosecond)

        return new_mdstate

    def clean_up(self):
        pass
        # shutil.rmtree(self.opt['out_directory'], ignore_errors=True)

        return
