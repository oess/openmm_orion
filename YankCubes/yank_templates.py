# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
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


yank_solvation_template = """\
---
options:
  verbose: {verbose}
  minimize: {minimize}
  output_dir: {output_directory}
  timestep: {timestep:f}*femtoseconds
  nsteps_per_iteration: {nsteps_per_iteration:d}
  number_of_iterations: {number_iterations:d}
  temperature: {temperature:f}*kelvin
  pressure: {pressure:f}*atmosphere
  anisotropic_dispersion_cutoff: auto
  resume_simulation: {resume_sim}
  resume_setup: {resume_sim}
  hydrogen_mass: {hydrogen_mass:f}*amu

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
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
        0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.35, 0.30, 
        0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00]

experiments:
  system: solvation-system
  protocol: solvation-protocol
"""

yank_binding_template = """\
---
options:
  verbose: {verbose}
  minimize: {minimize}
  output_dir: {output_directory}
  timestep: {timestep:f}*femtoseconds
  nsteps_per_iteration: {nsteps_per_iteration:d}
  number_of_iterations: {number_iterations:d}
  temperature: {temperature:f}*kelvin
  pressure: {pressure:f}*atmosphere
  anisotropic_dispersion_cutoff: auto
  resume_simulation: {resume_sim}
  resume_setup: {resume_sim}
  hydrogen_mass: {hydrogen_mass:f}*amu
  
systems:
  solvation-system:
    phase1_path: [{complex_pdb_fn}, {complex_xml_fn}]
    phase2_path: [{solvent_pdb_fn}, {solvent_xml_fn}]
    ligand_dsl: resname {ligand_resname}
    solvent_dsl: resname {solvent_dsl}

protocols:
  solvation-protocol:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.90, 0.78, 0.64, 0.51, 
        0.35, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 
        0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        lambda_restraints:     [0.00, 0.025, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00, 
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00] 
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.90, 0.78, 0.64, 0.51, 0.35, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
        0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.35, 0.30, 
        0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
experiments:
  system: solvation-system
  protocol: solvation-protocol
  restraint:
    type: {restraints}
"""