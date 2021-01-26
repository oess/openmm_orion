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

gromacs_nes = """
; RUN CONTROL PARAMETERS
integrator               = sd
dt                       = 0.002		; md timestep
nsteps                   = {nsteps:d}		; number of md steps
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 1000
nstvout                  = 0
nstfout                  = 0
nstlog                   = 10000
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
Pcoupl                   = Parrinello-Rahman
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
lincs-warnangle          = 30
morse                    = no

; FREE ENERGY CONTROL
free-energy              = yes
init-lambda              = 0
delta-lambda             = {dlambda:.2e}
sc-alpha                 = 0.3
sc-sigma                 = 0.25
sc-power                 = 1
sc-coul = yes
nstdhdl                  = 1
nstcalcenergy            = 1


; Non-equilibrium MD stuff
acc-grps                 = 
accelerate               = 
freezegrps               = ;FREEZE
freezedim                = ;Y Y Y
cos-acceleration         = 0
"""
