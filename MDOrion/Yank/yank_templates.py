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


# Total max running time per cube in hours
max_cube_running_time = 10.0

resources = {'k80':
                 {'w29': {'slope': 1.12e-6, 'intercept': 0.002},
                  'wsams': {'slope': 8.10e-8, 'intercept': 0.005}
                  }
             }  # intercept in hrs

yank_solvation_template = """
---
options:
  verbose: {verbose}
  minimize: {minimize}
  output_dir: {output_directory}
  temperature: {temperature:f}*kelvin
  pressure: {pressure:f}*atmosphere
  anisotropic_dispersion_cutoff: auto
  resume_simulation: {resume_sim}
  resume_setup: {resume_sim}
  hydrogen_mass: {hydrogen_mass:f}*amu
  processes_per_experiment: 1
  alchemical_pme_treatment: {alchemical_pme_treatment}
  checkpoint_interval: {checkpoint_interval}

mcmc_moves:
  langevin:
    type: LangevinSplittingDynamicsMove
    timestep: {timestep:f}*femtoseconds
    splitting: 'V R R R O R R R V'
    n_steps: {nsteps_per_iteration:d}

samplers:
  repex:
    type: ReplicaExchangeSampler
    mcmc_moves: langevin
    number_of_iterations: {number_iterations:d}
    online_analysis_interval: null

solvents:
  solvent:
    nonbonded_method: PME
    nonbonded_cutoff: 9*angstroms
    clearance: 8*angstroms
  vacuum:
    nonbonded_method: NoCutoff

systems:
  solvation-system:
    phase1_path: [{solvated_pdb_fn}, {solvated_xml_fn}]
    phase2_path: [{solute_pdb_fn}, {solute_xml_fn}]
    solvent1: solvent
    solvent2: vacuum
    solvent_dsl: resname {solvent_dsl}

protocols:
  solvation-protocol:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
                                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 
                                0.50, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00]

experiments:
  system: solvation-system
  protocol: solvation-protocol
  sampler: repex
"""


yank_binding_template = """
---
options:
  verbose: {verbose}
  minimize: {minimize}
  output_dir: {output_directory}
  temperature: {temperature:f}*kelvin
  pressure: {pressure:f}*atmosphere
  anisotropic_dispersion_cutoff: auto
  resume_simulation: {resume_sim}
  resume_setup: {resume_sim}
  hydrogen_mass: {hydrogen_mass:f}*amu
  processes_per_experiment: 1
  alchemical_pme_treatment: {alchemical_pme_treatment}
  checkpoint_interval: {checkpoint_interval}

mcmc_moves:
  langevin:
    type: LangevinSplittingDynamicsMove
    timestep: {timestep:f}*femtoseconds
    splitting: 'V R R R O R R R V'
    n_steps: {nsteps_per_iteration:d}

systems:
  system:
    phase1_path: [{complex_pdb_fn}, {complex_xml_fn}]
    phase2_path: [{solvent_pdb_fn}, {solvent_xml_fn}]
    ligand_dsl: resname {ligand_resname}
    solvent_dsl: resname {solvent_dsl}

samplers:
  repex:
    type: ReplicaExchangeSampler
    mcmc_moves: langevin
    number_of_iterations: {number_iterations:d}
    online_analysis_interval: null

  sams:
    type: SAMSSampler
    mcmc_moves: langevin
    state_update_scheme: global-jump
    flatness_threshold: 10.0
    number_of_iterations: {number_iterations:d}
    gamma0: 10.0
    online_analysis_interval: null

protocols:
  auto_protocol:
    solvent1:
      alchemical_path: auto
    solvent2:
      alchemical_path: auto

  windows_sams:
    solvent1:
      alchemical_path:
        lambda_restraints:     [0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00]
        lambda_electrostatics: [1.00, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.85, 0.80,
                0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.000, 0.00, 0.000, 0.00, 0.000, 0.00, 0.000, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.99,
                0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82,
                0.81, 0.80, 0.79, 0.78, 0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.70, 0.69, 0.68, 0.67, 0.66, 0.65,
                0.64, 0.63, 0.62, 0.61, 0.60, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52, 0.51, 0.50, 0.49, 0.48,
                0.47, 0.46, 0.45, 0.44, 0.43, 0.42, 0.41, 0.40, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.32, 0.31,
                0.30, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14,
                0.13, 0.12, 0.11, 0.10, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.05, 0.04, 0.03, 0.02,
                0.01, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.85, 0.80,
                0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.000, 0.00, 0.000, 0.00, 0.000, 0.00, 0.000, 0.00, 0.00, 0.00, 0.00, 0.00,
                0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.99,
                0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82,
                0.81, 0.80, 0.79, 0.78, 0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.70, 0.69, 0.68, 0.67, 0.66, 0.65,
                0.64, 0.63, 0.62, 0.61, 0.60, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52, 0.51, 0.50, 0.49, 0.48,
                0.47, 0.46, 0.45, 0.44, 0.43, 0.42, 0.41, 0.40, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.32, 0.31,
                0.30, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14,
                0.13, 0.12, 0.11, 0.10, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.05, 0.04, 0.03, 0.02,
                0.01, 0.00]
                
  windows_29:
    solvent1:
      alchemical_path:
        lambda_restraints:     [0.00, 0.025, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00,
                               1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                               1.00]
        lambda_electrostatics: [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.90, 0.78, 0.64, 0.51,
                                0.35, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                                0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05,
                                0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
                                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 
                                0.50, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]    
                          
harmonic:
  sampler: {sampler}
  system: system
  protocol: {protocol}
  restraint:
    type: Harmonic

boresch:
  sampler: {sampler}
  system: system
  protocol: {protocol}
  restraint:
    type: PeriodicTorsionBoresch

experiments: [{restraints}]
"""