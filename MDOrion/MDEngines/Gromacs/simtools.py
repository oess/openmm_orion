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

from simtk import unit

from simtk.openmm import app

import simtk

from MDOrion.MDEngines.utils import (MDSimulations,
                                     md_keys_converter)

import numpy as np

from MDOrion.MDEngines.Gromacs.gromacs_templates import (gromacs_minimization,
                                                         gromacs_nvt_npt)
import subprocess

import parmed

import copy

from oeommtools import utils as oeommutils

import tarfile

import os

import io

from openeye import oechem

from MDOrion.Standards import MDEngines, MDFileNames


class GromacsSimulations(MDSimulations):

    def __init__(self, mdstate,parmed_structure, opt):
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
        #
        # for c in topology.chains():
        #     for r in c.residues():
        #         for a in r.atoms():
        #             if r.name + a.name in atom_types_dic:
        #                 a.id = atom_types_dic[r.name + a.name]
        #             else:
        #                 a.id = 'O' + str(count_id)
        #                 count_id += 1
        #                 atom_types_dic[r.name + a.name] = a.id

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

        self.stepLen = 0.002 * unit.picoseconds

        opt['timestep'] = self.stepLen

        cutoff = opt['nonbondedCutoff'] * unit.angstroms

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
            cutoff = 0.0
            nslist = 0

        if opt['SimType'] == 'min':

            if opt['steps'] == 0:
                max_minimization_steps = 100000
            else:
                max_minimization_steps = opt['steps']

            mdp_template = gromacs_minimization.format(
                nsteps=max_minimization_steps,
                nslist=1,
                cutoff=cutoff.in_units_of(unit.nanometer) / unit.nanometer,
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

            if opt['trajectory_interval']:
                trajectory_steps = int(round(opt['trajectory_interval'] / (
                        opt['timestep'].in_units_of(unit.nanoseconds) / unit.nanoseconds)))
            else:
                trajectory_steps = 0

            # Convert simulation time in steps
            opt['steps'] = int(round(opt['time'] / (self.stepLen.in_units_of(unit.nanoseconds) / unit.nanoseconds)))

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
        opt['grm_gro_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".gro")
        opt['grm_pdb_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".pdb")
        opt['grm_tpr_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".tpr")
        opt['mdp_fn'] = os.path.join(opt['out_directory'], opt['outfname']+".mdp")
        opt['grm_def_fn'] = os.path.join(opt['out_directory'], opt['outfname']+"_run")
        opt['grm_log_fn'] = opt['grm_def_fn'] + '.log'
        opt['grm_trj_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".trr")
        opt['grm_trj_comp_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + ".xtc")

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

            protein_list = sorted(oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask='protein'))

            if protein_list:
                protein_range = [min(protein_list), max(protein_list)]
                last_idx = protein_range[1]
            else:
                protein_range = None

            ligand_list = sorted(oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask='ligand'))

            if ligand_list:
                ligand_range = [min(ligand_list), max(ligand_list)]
                last_idx =ligand_range[1]
            else:
                ligand_range = None

            water_list = sorted(oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask='water'))

            if water_list:
                water_range = [min(water_list), max(water_list)]
                last_idx = water_range[1]

            else:
                water_range = None

            exc_range = [last_idx + 1]

            if res_atom_list:
                apply_restraints = True

                restraints_dic = {}

                for atmol in opt['molecule'].GetAtoms():
                    at_idx = atmol.GetIdx()
                    if at_idx in res_atom_list:
                        res = oechem.OEAtomGetResidue(atmol)
                        if oechem.OEIsStandardProteinResidue(res):
                            if 'system1' in restraints_dic:
                                restraints_dic['system1'].append(at_idx)
                            else:
                                restraints_dic['system1'] = []
                                restraints_dic['system1'].append(at_idx)
                        else:
                            res_name = res.GetName()
                            if res_name in restraints_dic:
                                restraints_dic[res_name].append(at_idx)
                            else:
                                restraints_dic[res_name] = []
                                restraints_dic[res_name].append(at_idx)

        # Apply restraints
        if apply_restraints:

            # Generate topology files
            new_system_structure.save(opt['grm_top_fn'], overwrite=True)

            opt['Logger'].info("[{}] RESTRAINT mask applied to: {}"
                               "\tRestraint weight: {}".format(opt['CubeTitle'],
                                                               opt['restraints'],
                                                               opt[
                                                                   'restraintWt'] * unit.kilocalories_per_mole / unit.angstroms ** 2))

            # Select atom to restraint
            opt['Logger'].info("[{}] Number of restraint atoms: {}".format(opt['CubeTitle'], len(res_atom_list)))

            # Restrained atom index file to be appended to the index file
            opt['grm_res_idx_fn'] = os.path.join(opt['out_directory'], opt['outfname'] + '_res.idx')

            digits = len(str(abs(res_atom_list[-1])))

            # Gromacs number of column in the restraint file
            chunk = 15

            count = 0
            with open(opt['grm_res_idx_fn'], 'w') as f:

                for key in restraints_dic.keys():
                    f.write("[ {} ]\n".format(key))
                    for i in range(0, len(restraints_dic[key])):

                        sel = restraints_dic[key][i]

                        if protein_range is not None and protein_range[0] <= sel <= protein_range[1]:
                            index_rs = sel - protein_range[0]
                        elif ligand_range is not None and ligand_range[0] <= sel <= ligand_range[1]:
                            index_rs = sel - ligand_range[0]
                        elif water_range is not None and water_range[0] <= sel <= water_range[1]:
                            index_rs = sel - water_range[0]
                        else:
                            index_rs = sel - exc_range[0]

                        f.write("{:>{digits}}".format(str(index_rs + 1), digits=digits))
                        if count == chunk - 1 or i == len(res_atom_list) - 1:
                            count = 0
                            f.write("\n")
                        else:
                            f.write(" ")
                            count += 1
                    f.write('\n')

            # Restrained System index file
            opt['grm_system_fn'] = os.path.join(opt['out_directory'], opt['outfname']+'_sys')

            p = subprocess.Popen(['gmx',
                                  'make_ndx',
                                  '-f', opt['grm_gro_fn'],
                                  '-o', opt['grm_system_fn']],
                                 stdin=subprocess.PIPE)

            # Select the System
            p.communicate(b'0\nq\n')

            append_fns = [opt['grm_system_fn']+'.ndx', opt['grm_res_idx_fn']]

            # Index file name
            opt['grm_ndx_fn'] = os.path.join(opt['out_directory'], opt['outfname']+'.ndx')

            with open(opt['grm_ndx_fn'], 'w') as outfile:
                for fn in append_fns:
                    with open(fn) as infile:
                        outfile.write(infile.read())

            with open(opt['grm_ndx_fn'], 'r') as f:
                lines_res = f.readlines()

            count_systems = 0
            for line in lines_res:
                if '[' in line:
                    count_systems += 1

            # Restrains position file name
            opt['grm_itp_fn'] = os.path.join(opt['out_directory'], opt['outfname'])

            indx_selection = count_systems

            # Restraints weight
            res_wgt = opt['restraintWt'] * unit.kilocalories_per_mole / (unit.angstroms ** 2)

            # grm unit
            res_wgt_grm = str(int(res_wgt.in_units_of(unit.kilojoule_per_mole / (unit.nanometer ** 2)) / (
                        unit.kilojoule_per_mole / (unit.nanometer ** 2)))).encode()

            for key, idx in zip(restraints_dic.keys(), range(indx_selection - len(restraints_dic), indx_selection)):

                p = subprocess.Popen(['gmx', 'genrestr',
                                      '-f', opt['grm_gro_fn'],
                                      '-n', opt['grm_ndx_fn'],
                                      '-o', opt['grm_itp_fn']+'_'+key+'.itp',
                                      '-fc', res_wgt_grm, res_wgt_grm, res_wgt_grm],
                                     stdin=subprocess.PIPE)

                p.communicate(str(idx).encode())

            with open(opt['grm_top_fn'], 'r') as f:
                lines = f.readlines()

            # itp file names
            itp_include_fns = {key: opt['grm_itp_fn']+'_'+key+'.itp' for key in restraints_dic.keys()}

            for key in itp_include_fns.keys():
                for idx in range(0, len(lines)):
                    if lines[idx] == "[ moleculetype ]\n" and lines[idx + 2].startswith(key):
                        for count in range(idx + 2, len(lines)):
                            if lines[count] == "[ moleculetype ]\n":
                                include = """                
; Include Position restraint file
#ifdef POSRES
#include "{}"
#endif\n\n
""".format(itp_include_fns[key])
                                lines.insert(count - 1, include)

                                break
                        break

            with open(opt['grm_top_fn'], 'w') as f:
                for line in lines:
                    f.write(line)

        # Generate Gromacs .tpr file
        if apply_restraints:
            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', opt['mdp_fn'],
                                   '-c', opt['grm_gro_fn'],
                                   '-r', opt['grm_gro_fn'],
                                   '-p', opt['grm_top_fn'],
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
        subprocess.check_call(['gmx',
                               'mdrun',
                               '-v',
                               '-s', self.opt['grm_tpr_fn'],
                               '-deffnm', self.opt['grm_def_fn'],
                               '-o', self.opt['grm_trj_fn']
                               ])

        if self.opt['SimType'] in ['nvt', 'npt']:

            if self.opt['reporter_interval']:

                with(io.open(self.opt['grm_log_fn'], 'r', encoding='utf8', errors='ignore')) as fr:
                    log_string = fr.read()

                self.opt['str_logger'] += '\n'+log_string

            # Save trajectory files
            if self.opt['trajectory_interval']:

                # Generate whole system trajectory
                p = subprocess.Popen(['gmx',
                                      'trjconv',
                                      '-f', self.opt['grm_trj_fn'],
                                      '-s', self.opt['grm_tpr_fn'],
                                      '-o', self.opt['grm_trj_comp_fn'],
                                      '-pbc', b'whole'],
                                     stdin=subprocess.PIPE)

                # Select the entire System
                p.communicate(b'0')

                # Tar the files dir with its content:
                tar_fn = self.opt['trj_fn']

                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(self.opt['grm_gro_fn'], arcname=os.path.basename(self.opt['grm_gro_fn']))
                    archive.add(self.opt['grm_pdb_fn'], arcname=os.path.basename(self.opt['grm_pdb_fn']))
                    archive.add(self.opt['grm_top_fn'], arcname=os.path.basename(self.opt['grm_top_fn']))
                    archive.add(self.opt['grm_trj_comp_fn'], arcname=os.path.basename(self.opt['grm_trj_comp_fn']))
                    archive.add(self.opt['grm_log_fn'], arcname=os.path.basename(self.opt['grm_log_fn']))

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
