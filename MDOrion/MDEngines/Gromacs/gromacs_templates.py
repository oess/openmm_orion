gromacs_minimization = """
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = {nsteps:d}         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme   = Verlet    ; Buffered neighbor searching
nstlist		    = {nslist:d} ;  Frequency to update the neighbor list and long range forces
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = {cutoff:f}       ; Short-range electrostatic cut-off
rvdw            = {cutoff:f}       ; Short-range Van der Waals cut-off
pbc             = {pbc}     ; Periodic Boundary Conditions 
"""


gromacs_nvt_npt = """
define		= -DPOSRES	; position restrain

; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= {nsteps:d}		; number od md steps
dt		    = {timestep:f}		; md timestep
comm-mode   = Linear ;  mode for center of mass motion removal
nstcomm     = 100 ;number of steps for center of mass motion removal

; Output control
nstxout		= {trajectory_steps:d}		; save coordinates
nstvout		= {trajectory_steps:d}		; save velocities
nstenergy	= {trajectory_steps:d}		; save energies
nstcalcenergy = {trajectory_steps:d}    ; calculate energy every provided steps
nstlog		= {reporter_steps:d}		; update log file

; Bond parameters
continuation	        = no		; Restarting
constraint_algorithm    = lincs	    ; holonomic constraints
constraints	            = {constraints}	; constraint type
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy

; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = {nslist:d}	    ; largely irrelevant with Verlet scheme
rcoulomb	    = {cutoff:f}		; short-range electrostatic cutoff (in nm)
rvdw		    = {cutoff:f}		; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT

; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= System	;
tau_t		= 0.1	    ; time constant, in ps
ref_t		= {temperature:f} ; reference temperature in K

; Pressure coupling is on
pcoupl		        = {pcoupl}	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 0.5		            ; time constant, in ps
ref_p		        = {pressure:f}		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com

; Periodic boundary conditions
pbc		= {pbc}		; Periodic Boundary Conditions

; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme

; Velocity generation
gen_vel		= {gen_vel}		; Velocity generation
gen_temp    = {temperature:f} ; reference temperature in K
"""

gromacs_pos_restraints = "{:<{digits}}\t{:<}  {:<5.2f}  {:<5.2f}  {:<5.2f}\n"

# LANGEVIN DYNAMICS
# gromacs_nvt_npt = """
# define		= -DPOSRES	; position restrain
#
# ; Run parameters
# integrator	= sd		; leap-frog integrator
# nsteps		= {nsteps:d}		; number of md steps
# dt		    = {timestep:f}		; md timestep
#
# ; Output control
# nstxout		= {trajectory_steps:d}		; save coordinates
# nstvout		= {trajectory_steps:d}		; save velocities
# nstenergy	= {trajectory_steps:d}		; save energies
# nstcalcenergy = {trajectory_steps:d}    ; calculate energy every provided steps
# nstlog		= {reporter_steps:d}		; update log file
#
# ; Bond parameters
# continuation	        = no		; Restarting
# constraint_algorithm    = lincs	    ; holonomic constraints
# constraints	            = {constraints}	; constraint type
# lincs_iter	            = 1		    ; accuracy of LINCS
# lincs_order	            = 4		    ; also related to accuracy
#
# ; Neighborsearching
# cutoff-scheme   = Verlet
# ns_type		    = grid		; search neighboring grid cells
# nstlist		    = {nslist:d}	    ; largely irrelevant with Verlet scheme
# rcoulomb	    = {cutoff:f}		; short-range electrostatic cutoff (in nm)
# rvdw		    = {cutoff:f}		; short-range van der Waals cutoff (in nm)
#
# ; Electrostatics
# coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
# pme_order	    = 4		    ; cubic interpolation
# fourierspacing	= 0.16		; grid spacing for FFT
#
# ; Temperature coupling is on
# tcoupl		= no
# tc-grps		= System	;
# tau_t		= 0.1	    ; time constant, in ps
# ref_t		= {temperature:f} ; reference temperature in K
#
# ; Pressure coupling is on
# pcoupl		        = {pcoupl}	    ; Pressure coupling on in NPT
# pcoupltype	        = isotropic	            ; uniform scaling of box vectors
# tau_p		        = 0.5		            ; time constant, in ps
# ref_p		        = {pressure:f}		            ; reference pressure, in bar
# compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
# refcoord_scaling        = com
#
# ; Periodic boundary conditions
# pbc		= {pbc}		; Periodic Boundary Conditions
#
# ; Dispersion correction
# DispCorr	= EnerPres	; account for cut-off vdW scheme
#
# ; Velocity generation
# gen_vel		= {gen_vel}		; Velocity generation
# gen_temp    = {temperature:f} ; reference temperature in K
# """