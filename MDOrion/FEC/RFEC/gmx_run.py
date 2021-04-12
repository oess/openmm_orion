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


import subprocess

from subprocess import STDOUT, PIPE, Popen, DEVNULL

from simtk import unit

import os

from tempfile import TemporaryDirectory


gromacs_min_nes_nvt_npt = """
; RUN CONTROL PARAMETERS

define		= {restraints}	; position restrain

integrator               = {integrator}
dt                       = 0.002		; md timestep
nsteps                   = {nsteps:d}		; number of md steps
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; ENERGY MINIMIZATION OPTIONS
; Force tolerance and initial step-size
emtol                    = 100
emstep                   = 0.01
; Max number of iterations in relax_shells
niter                    = 0
; Step size (1/ps^2) for minimization of flexible constraints
fcstep                   = 0
; Frequency of steepest descents steps when doing CG
nstcgsteep               = 1000
nbfgscorr                = 10


; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 0
nstenergy                = 0

; NEIGHBORSEARCHING PARAMETERS
cutoff-scheme = verlet
verlet-buffer-tolerance = 0.00001
nstlist                  = 10
ns-type                  = Grid
pbc                      = xyz     
rlist                    = 1.2


; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = {cutoff:.2f}
epsilon-r                = 1
vdw-type                 = switch      
rvdw-switch              = {rvdwswitch:.2f}
rvdw                     = {cutoff:.2f}
DispCorr                 = EnerPres
table-extension          = 1
fourierspacing           = 0.12
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
pme_order                = 4
ewald_rtol               = 1e-05
ewald_geometry           = 3d
epsilon_surface          = 0


; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling  
tcoupl                   = v-rescale
; Groups to couple separately
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau-t                    = 2.0
ref-t                    = {temperature:f} ; reference temperature in K

; Pressure coupling     
Pcoupl                   = {barostat}
Pcoupltype               = Isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau-p                    = 5
compressibility          = 4.6E-5
ref-p                    = {pressure:f} ; reference pressure, in bar

; GENERATE VELOCITIES FOR STARTUP RUN
gen-vel                  = {gen_vel}
continuation	         = {continue_sim}

; OPTIONS FOR BONDS    
constraints              = {lincs_type}
constraint-algorithm     = lincs
lincs-order              = 4
lincs-iter               = 2
lincs-warnangle          = 50
morse                    = no

; FREE ENERGY CONTROL
free-energy              = yes
init-lambda              = 0
delta-lambda             = {dlambda:.2e}
sc-alpha                 = 0.3
sc-sigma                 = 0.25
sc-power                 = 1
sc-coul = yes
nstdhdl                  = {nstdhdl:d}
nstcalcenergy            = 1

domain-decomposition     = no

; Non-equilibrium MD stuff
acc-grps                 = 
accelerate               = 
freezegrps               = ;FREEZE
freezedim                = ;Y Y Y
cos-acceleration         = 0
"""


def check_gmx_grompp(gro, top, sim_type=None, verbose=False):
    try:
        with TemporaryDirectory() as outdir:

            gro_fn = os.path.join(outdir, "gmx_gro.gro")
            top_fn = os.path.join(outdir, "gmx_top.top")

            with open(gro_fn, 'w') as f:
                f.write(gro)
            with open(top_fn, 'w') as f:
                f.write(top)

            gmx_gro_fn = os.path.join(outdir, gro_fn)
            gmx_top_fn = os.path.join(outdir, top_fn)
            gmx_tpr_fn = os.path.join(outdir, "gmx_tpr.tpr")
            gmx_mdp_fn = os.path.join(outdir, "gmx_mdp.mdp")

            gmx_fe_template = gromacs_min_nes_nvt_npt.format(restraints='',
                                                             integrator='steep',
                                                             nsteps=300,
                                                             temperature=300,
                                                             barostat='Parrinello-Rahman',
                                                             pressure=1.1,
                                                             gen_vel='no',
                                                             continue_sim='no',
                                                             lincs_type='all-bonds',
                                                             cutoff=1.0,
                                                             rvdwswitch=0.9,
                                                             dlambda=0,
                                                             nstdhdl=0)

            with open(gmx_mdp_fn, 'w') as f:
                f.write(gmx_fe_template)

            if verbose:

                subprocess.check_call(['gmx',
                                       'grompp',
                                       '-f', gmx_mdp_fn,
                                       '-c', gmx_gro_fn,
                                       '-p', gmx_top_fn,
                                       '-o', gmx_tpr_fn,
                                       '-maxwarn', '4'
                                       ])

                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', gmx_tpr_fn])

            else:

                subprocess.check_call(['gmx',
                                       'grompp',
                                       '-f', gmx_mdp_fn,
                                       '-c', gmx_gro_fn,
                                       '-p', gmx_top_fn,
                                       '-o', gmx_tpr_fn,
                                       '-maxwarn', '4'
                                       ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT)

                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', gmx_tpr_fn,
                                       ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT, timeout=300)

    except:
        raise ValueError("Cannot Assemble the Gromacs .tpr file for the sim type: {}".format(sim_type))


def _run_gmx(mdp_fn, gro_fn, top_fn, tpr_fn, deffnm_fn, opt, cpti_fn=None, cpto_fn=None):

    if opt['verbose']:

        # Assemble the Gromacs system to run
        if not opt['restraints']:
            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', mdp_fn,
                                   '-c', gro_fn,
                                   '-p', top_fn,
                                   '-o', tpr_fn,
                                   '-maxwarn', '4'
                                   ])
        else:
            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', mdp_fn,
                                   '-c', gro_fn,
                                   '-r', gro_fn,
                                   '-p', top_fn,
                                   '-o', tpr_fn,
                                   '-maxwarn', '4'
                                   ])

        # Run Gromacs
        if cpti_fn is None and cpto_fn is None:
            subprocess.check_call(['gmx',
                                   'mdrun',
                                   '-v',
                                   '-s', tpr_fn,
                                   '-deffnm', deffnm_fn
                                   ])
        else:

            if cpti_fn is None and cpto_fn is not None:

                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', tpr_fn,
                                       '-cpo', cpto_fn,
                                       '-deffnm', deffnm_fn
                                       ])
            else:
                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', tpr_fn,
                                       '-noappend',
                                       '-cpi', cpti_fn,
                                       '-cpo', cpto_fn,
                                       '-deffnm', deffnm_fn
                                       ])

    else:
        # Assemble the Gromacs system to run
        if not opt['restraints']:
            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', mdp_fn,
                                   '-c', gro_fn,
                                   '-p', top_fn,
                                   '-o', tpr_fn,
                                   '-maxwarn', '4'
                                   ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT,
                                  timeout=opt['gmx_process_timeout'])
        else:
            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', mdp_fn,
                                   '-c', gro_fn,
                                   '-r', gro_fn,
                                   '-p', top_fn,
                                   '-o', tpr_fn,
                                   '-maxwarn', '4'
                                   ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT,
                                  timeout=opt['gmx_process_timeout'])

        # Run Gromacs
        if cpti_fn is None and cpto_fn is None:
            subprocess.check_call(['gmx',
                                   'mdrun',
                                   '-v',
                                   '-s', tpr_fn,
                                   '-deffnm', deffnm_fn
                                   ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT,
                                  timeout=opt['gmx_process_timeout'])
        else:

            if cpti_fn is None and cpto_fn is not None:
                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', tpr_fn,
                                       '-cpo', cpto_fn,
                                       '-deffnm', deffnm_fn
                                       ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT,
                                      timeout=opt['gmx_process_timeout'])
            else:
                subprocess.check_call(['gmx',
                                       'mdrun',
                                       '-v',
                                       '-s', tpr_fn,
                                       '-noappend',
                                       '-cpi', cpti_fn,
                                       '-cpo', cpto_fn,
                                       '-deffnm', deffnm_fn
                                       ], stdin=PIPE, stdout=DEVNULL, stderr=STDOUT,
                                      timeout=opt['gmx_process_timeout'])


def gmx_run(gmx_gro, gmx_top, opt):

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

    out_dir = opt['out_directory']

    # TODO DEBUGGING
    # out_dir = "./"

    # FAST MINIMIZATION
    nsteps = 2000
    opt['restraints'] = True

    gmx_gro_fn = os.path.join(out_dir, "gmx_gro.gro")
    gmx_top_fn = os.path.join(out_dir, "gmx_top.top")
    gmx_min_mdp_fn = os.path.join(out_dir, "gmx_min.mdp")
    gmx_min_tpr_fn = os.path.join(out_dir, "gmx_min.tpr")
    gmx_min_deffnm_fn = os.path.join(out_dir, "gmx_run_min_" + str(opt['frame_count']))
    gmx_min_confout_gro = gmx_min_deffnm_fn + '.gro'

    with open(gmx_gro_fn, 'w') as f:
        f.write(gmx_gro)

    with open(gmx_top_fn, 'w') as f:
        f.write(gmx_top)

    gmx_fe_template = gromacs_min_nes_nvt_npt.format(restraints='-DRESTRAINTS',
                                                     integrator='steep',
                                                     nsteps=nsteps,
                                                     temperature=opt['temperature'],
                                                     barostat='no',
                                                     pressure=pressure.value_in_unit(unit.bar),
                                                     gen_vel='no',
                                                     continue_sim='no',
                                                     lincs_type=opt['lincs_type'],
                                                     cutoff=cutoff_distance.value_in_unit(unit.nanometer),
                                                     rvdwswitch=rvdw_switch.value_in_unit(unit.nanometer),
                                                     dlambda=0,
                                                     nstdhdl=0)
    with open(gmx_min_mdp_fn, 'w') as f:
        f.write(gmx_fe_template)

    _run_gmx(gmx_min_mdp_fn, gmx_gro_fn, gmx_top_fn, gmx_min_tpr_fn, gmx_min_deffnm_fn, opt,
             cpti_fn=None, cpto_fn=None)

    # Restrained NVT
    nsteps = 5000
    opt['restraints'] = True

    gmx_eq_nvt_mdp_fn = os.path.join(out_dir, "gmx_eq_nvt.mdp")
    gmx_eq_nvt_tpr_fn = os.path.join(out_dir, "gmx_eq_nvt.tpr")
    gmx_nvt_deffnm_fn = os.path.join(out_dir, "gmx_run_nvt_" + str(opt['frame_count']))
    gmx_nvt_cpto_fn = gmx_nvt_deffnm_fn + '.cpt'
    gmx_nvt_confout_gro = gmx_nvt_deffnm_fn + '.gro'

    # gmx_trj_fn = gmx_nvt_deffnm_fn + '.trr'

    gmx_fe_template = gromacs_min_nes_nvt_npt.format(restraints='-DRESTRAINTS',
                                                     integrator='md',
                                                     nsteps=nsteps,
                                                     temperature=opt['temperature'],
                                                     barostat='no',
                                                     pressure=pressure.value_in_unit(unit.bar),
                                                     gen_vel='no',
                                                     continue_sim='yes',
                                                     lincs_type=opt['lincs_type'],
                                                     cutoff=cutoff_distance.value_in_unit(unit.nanometer),
                                                     rvdwswitch=rvdw_switch.value_in_unit(unit.nanometer),
                                                     dlambda=0,
                                                     nstdhdl=0)

    with open(gmx_eq_nvt_mdp_fn, 'w') as f:
        f.write(gmx_fe_template)

    _run_gmx(gmx_eq_nvt_mdp_fn, gmx_min_confout_gro, gmx_top_fn, gmx_eq_nvt_tpr_fn, gmx_nvt_deffnm_fn, opt,
             cpti_fn=None, cpto_fn=None)

    # Restrained NPT
    nsteps = 5000
    opt['restraints'] = True

    gmx_eq_npt0_mdp_fn = os.path.join(out_dir, "gmx_eq_npt0.mdp")
    gmx_eq_npt0_tpr_fn = os.path.join(out_dir, "gmx_eq_npt0.tpr")
    gmx_npt0_deffnm_fn = os.path.join(out_dir, "gmx_run_npt0" + '_' + str(opt['frame_count']))
    gmx_npt0_cpto_fn = gmx_npt0_deffnm_fn + '.cpt'
    gmx_npt0_confout_gro = gmx_npt0_deffnm_fn + '.gro'

    gmx_fe_template = gromacs_min_nes_nvt_npt.format(restraints='-DRESTRAINTS',
                                                     integrator='md',
                                                     nsteps=nsteps,
                                                     temperature=opt['temperature'],
                                                     barostat='Parrinello-Rahman',
                                                     pressure=pressure.value_in_unit(unit.bar),
                                                     gen_vel='no',
                                                     continue_sim='yes',
                                                     lincs_type=opt['lincs_type'],
                                                     cutoff=cutoff_distance.value_in_unit(unit.nanometer),
                                                     rvdwswitch=rvdw_switch.value_in_unit(unit.nanometer),
                                                     dlambda=0,
                                                     nstdhdl=0)

    with open(gmx_eq_npt0_mdp_fn, 'w') as f:
        f.write(gmx_fe_template)

    _run_gmx(gmx_eq_npt0_mdp_fn, gmx_nvt_confout_gro, gmx_top_fn, gmx_eq_npt0_tpr_fn, gmx_npt0_deffnm_fn, opt,
             cpti_fn=None, cpto_fn=gmx_npt0_cpto_fn)

    # NPT
    opt['restraints'] = False
    nsteps = 5000

    gmx_eq_npt1_mdp_fn = os.path.join(out_dir, "gmx_eq_npt1.mdp")
    gmx_eq_npt1_tpr_fn = os.path.join(out_dir, "gmx_eq_npt1.tpr")
    gmx_npt1_deffnm_fn = os.path.join(out_dir, "gmx_run_npt1" + '_' + str(opt['frame_count']))
    gmx_npt1_cpti_fn = gmx_npt0_cpto_fn
    gmx_npt1_cpto_fn = gmx_npt1_deffnm_fn + '.cpt'
    gmx_npt1_confout_gro = gmx_npt1_deffnm_fn + '.gro'

    gmx_fe_template = gromacs_min_nes_nvt_npt.format(restraints='',
                                                     integrator='md',
                                                     nsteps=nsteps,
                                                     temperature=opt['temperature'],
                                                     barostat='Parrinello-Rahman',
                                                     pressure=pressure.value_in_unit(unit.bar),
                                                     gen_vel='no',
                                                     continue_sim='yes',
                                                     lincs_type=opt['lincs_type'],
                                                     cutoff=cutoff_distance.value_in_unit(unit.nanometer),
                                                     rvdwswitch=rvdw_switch.value_in_unit(unit.nanometer),
                                                     dlambda=0,
                                                     nstdhdl=0)

    with open(gmx_eq_npt1_mdp_fn, 'w') as f:
        f.write(gmx_fe_template)

    _run_gmx(gmx_eq_npt1_mdp_fn, gmx_npt0_confout_gro, gmx_top_fn, gmx_eq_npt1_tpr_fn, gmx_npt1_deffnm_fn, opt,
             cpti_fn=None,
             cpto_fn=gmx_npt1_cpto_fn)

    # NES
    opt['restraints'] = False
    stepLen = 2.0 * unit.femtoseconds
    nsteps = int(round(opt['time'] / (stepLen.in_units_of(unit.nanoseconds) / unit.nanoseconds)))

    # Full decoupling lambda in [0,1]
    dlambda = 1.0 / nsteps

    gmx_ne_mdp_fn = os.path.join(out_dir, "gmx_ne.mdp")
    gmx_ne_tpr_fn = os.path.join(out_dir, "gmx_ne.tpr")
    gmx_ne_deffnm_fn = os.path.join(out_dir, "gmx_run_ne" + '_' + str(opt['frame_count']))
    gmx_ne_cpti_fn = gmx_npt1_cpto_fn
    gmx_ne_cpto_fn = gmx_ne_deffnm_fn + '.cpt'
    gmx_ne_confout_gro = gmx_ne_deffnm_fn + '.gro'

    gmx_fe_template = gromacs_min_nes_nvt_npt.format(restraints='',
                                                     integrator='md',
                                                     nsteps=nsteps,
                                                     temperature=opt['temperature'],
                                                     barostat='Parrinello-Rahman',
                                                     pressure=pressure.value_in_unit(unit.bar),
                                                     gen_vel='no',
                                                     continue_sim='yes',
                                                     lincs_type=opt['lincs_type'],
                                                     cutoff=cutoff_distance.value_in_unit(unit.nanometer),
                                                     rvdwswitch=rvdw_switch.value_in_unit(unit.nanometer),
                                                     dlambda=dlambda,
                                                     nstdhdl=1)

    with open(gmx_ne_mdp_fn, 'w') as f:
        f.write(gmx_fe_template)

    _run_gmx(gmx_ne_mdp_fn, gmx_npt1_confout_gro, gmx_top_fn, gmx_ne_tpr_fn, gmx_ne_deffnm_fn, opt,
             cpti_fn=None, cpto_fn=gmx_ne_cpto_fn)

    return gmx_ne_deffnm_fn
