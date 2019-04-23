#!/usr/bin/env python

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


from floe.api import WorkFloe

from MDOrion.MDEngines.cubes import (MDMinimizeCube,
                                     MDNvtCube,
                                     MDNptCube)

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.System.cubes import SolvationCube

from MDOrion.ProtPrep.cubes import ProteinSetting

from MDOrion.ForceField.cubes import ForceFieldCube

from MDOrion.LigPrep.cubes import (LigandChargeCube,
                                   LigandSetting)

from MDOrion.System.cubes import IDSettingCube

from MDOrion.Yank.cubes import (SyncBindingFECube,
                                YankBindingFECube,
                                YankProxyCube)

from cuberecord import (DatasetWriterCube,
                        DatasetReaderCube)

from MDOrion.TrjAnalysis.cubes import MDFloeReportCube

job = WorkFloe('Binding Affinity Self-Adjusted Mixture Sampling',
               title='Binding Affinity Self-Adjusted Mixture Sampling')

job.description = """
NOTE: this is an Alpha Test version. 
We are actively working on improving Yank in Orion

The Absolute Binding Affinity Free Energy protocol (ABFE) performs Binding Affinity calculations
on a set of provided ligands and related receptor by using YANK Self-Adjusted Mixture Sampling
( http://getyank.org/latest/ ). The ligands need to have coordinates and correct chemistry.
Each ligand can have multiple conformers, but each conformer will be treated as a different ligand
and prepared to run ABFE. The protein needs to be prepared at MD preparation standard. This includes
capping the protein, resolve missing atoms in protein residues and resolve missing protein loops.
The parametrization of some "known unknown" non standard protein residues is partially supported.
Ligands need to be already posed in the protein binding site. A complex (Bonded State) is formed,
solvated and parametrized accordingly to the selected force fields. In a similar fashion the Unbounded
state is also prepared. Minimization, Warm up (NVT) and Equilibration (NPT) stages are performed
an the Bonded and Unbounded states. In order to minimize Molecular Dynamics (MD) failures along
these stages, positional harmonic restraints are applied on the ligand and protein with different
force constants. At the end of the equilibration stages the ABFE calculations are run by YANK with
the selected parameters. Calculated Binding Affinities for each ligand are output with the related
floe reports.

Required Input Parameters:
-----------
ligands: Dataset of the prepared ligands
protein: Dataset of the prepared protein

Outputs:
--------
out : Dataset of the solvated systems with the calculated binding free energies and
floe reports
"""

job.classification = [['BindingFreeEnergy', 'Yank']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")
job.add_cube(iligs)

chargelig = LigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)
job.add_cube(chargelig)

ligset = LigandSetting("LigandSetting", title="Ligand Setting")
ligset.set_parameters(lig_res_name='LIG')
job.add_cube(ligset)

ligid = IDSettingCube("Ligand Ids")
job.add_cube(ligid)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = DatasetReaderCube("ProteinReader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                        description="Protein file name")
job.add_cube(iprot)

protset = ProteinSetting("ProteinSetting", title="Protein Setting")
protset.promote_parameter("protein_prefix", promoted_name="protein_prefix", default="PRT")
job.add_cube(protset)

# COMPLEX SETTING

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex", title="Complex Preparation")
complx.set_parameters(lig_res_name='LIG')
job.add_cube(complx)

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvateComplex = SolvationCube("HydrationComplex", title="Complex Hydration")
solvateComplex.promote_parameter('density', promoted_name='density', default=1.03,
                                 description="Solution density in g/ml")
solvateComplex.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                                 description='Salt concentration (Na+, Cl-) in millimolar')
solvateComplex.set_parameters(solvents='[H]O[H]')
solvateComplex.set_parameters(molar_fractions='1.0')
solvateComplex.set_parameters(close_solvent=True)
job.add_cube(solvateComplex)

# Complex Force Field Application
ffComplex = ForceFieldCube("ForceFieldComplex", title="Complex Parametrization")
ffComplex.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber99SBildn')
ffComplex.promote_parameter('ligand_forcefield', promoted_name='ligand_forcefield', default='Gaff2')
ffComplex.promote_parameter('other_forcefield', promoted_name='other_forcefield', default='Gaff2')
ffComplex.set_parameters(lig_res_name='LIG')
job.add_cube(ffComplex)

# Add YANK Cube
yank_proxy = YankProxyCube("YankProxy", title="Yank Proxy")
yank_proxy.promote_parameter('iterations', promoted_name='iterations',
                             default=1000,
                             description="Total number of Yank iterations")
job.add_cube(yank_proxy)

# First Yank Cube used to build the UI interface
abfe = YankBindingFECube("YankABFE", title="Yank ABFE SAMS")
abfe.promote_parameter('iterations', promoted_name='iterations',
                       description="Total Number of Yank iterations for the entire floe. "
                                   "A Yank iteration is 500 MD steps")
abfe.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                       description='Temperature (Kelvin)')
abfe.promote_parameter('pressure', promoted_name='pressure', default=1.0,
                       description='Pressure (atm)')
abfe.promote_parameter('hmr', promoted_name='hmr', default=False,
                       description='Hydrogen Mass Repartitioning')
abfe.promote_parameter('restraints', promoted_name='restraints',
                        default='boresch',
                       description='Select the restraint types to apply to the ligand during the '
                                   'alchemical decoupling. Choices: harmonic, boresch')
abfe.set_parameters(lig_res_name='LIG')
abfe.promote_parameter('verbose', promoted_name='verbose', default=False, description="Yank verbose mode on/off")
# abfe.promote_parameter('user_yank_yaml_file', promoted_name='yaml', default=None)
abfe.set_parameters(sampler='sams')
abfe.promote_parameter('protocol_sams', promoted_name='protocol_sams', default='windows_sams')
job.add_cube(abfe)

# Minimization
minComplex = MDMinimizeCube('minComplex', title='Complex Minimization')
minComplex.promote_parameter('hmr', promoted_name='hmr')
minComplex.set_parameters(restraints="noh (ligand or protein)")
minComplex.set_parameters(restraintWt=5.0)
minComplex.set_parameters(steps=1000)
minComplex.set_parameters(center=True)
job.add_cube(minComplex)

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmupComplex = MDNvtCube('warmupComplex', title='Complex Warm Up')
warmupComplex.set_parameters(time=0.02)
warmupComplex.promote_parameter('temperature', promoted_name='temperature')
warmupComplex.promote_parameter('hmr', promoted_name='hmr')
warmupComplex.set_parameters(restraints="noh (ligand or protein)")
warmupComplex.set_parameters(restraintWt=2.0)
warmupComplex.set_parameters(trajectory_interval=0.0)
warmupComplex.set_parameters(reporter_interval=0.0)
warmupComplex.set_parameters(suffix='warmup_complex')
job.add_cube(warmupComplex)

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1Complex = MDNptCube('equil1Complex', title='Complex Equilibration I')
equil1Complex.set_parameters(time=0.02)
equil1Complex.promote_parameter('temperature', promoted_name='temperature')
equil1Complex.promote_parameter('pressure', promoted_name='pressure')
equil1Complex.promote_parameter('hmr', promoted_name='hmr')
equil1Complex.set_parameters(restraints="noh (ligand or protein)")
equil1Complex.set_parameters(restraintWt=2.0)
equil1Complex.set_parameters(trajectory_interval=0.0)
equil1Complex.set_parameters(reporter_interval=0.0)
equil1Complex.set_parameters(suffix='equil1')
job.add_cube(equil1Complex)

# NPT Equilibration stage 2
equil2Complex = MDNptCube('equil2Complex', title='Complex Equilibration II')
equil2Complex.set_parameters(time=0.02)
equil2Complex.promote_parameter('temperature', promoted_name='temperature')
equil2Complex.promote_parameter('pressure', promoted_name='pressure')
equil2Complex.promote_parameter('hmr', promoted_name='hmr')
equil2Complex.set_parameters(restraints="noh (ligand or protein)")
equil2Complex.set_parameters(restraintWt=0.5)
equil2Complex.set_parameters(trajectory_interval=0.0)
equil2Complex.set_parameters(reporter_interval=0.0)
equil2Complex.set_parameters(suffix='equil2')
job.add_cube(equil2Complex)

# NPT Equilibration stage 3
equil3Complex = MDNptCube('equil3Complex', title='Complex Equilibration III')
equil3Complex.set_parameters(time=0.02)
equil3Complex.promote_parameter('temperature', promoted_name='temperature')
equil3Complex.promote_parameter('pressure', promoted_name='pressure')
equil3Complex.promote_parameter('hmr', promoted_name='hmr')
equil3Complex.set_parameters(restraints="ca_protein or (noh ligand)")
equil3Complex.set_parameters(restraintWt=0.1)
equil3Complex.set_parameters(trajectory_interval=0.0)
equil3Complex.set_parameters(reporter_interval=0.0)
equil3Complex.set_parameters(suffix='equil3')
job.add_cube(equil3Complex)

# LIGAND SETTING

# Solvate Ligands
solvateLigand = SolvationCube("HydrationLigand", title="Unbound Ligand Hydration")
solvateLigand.set_parameters(solvents='[H]O[H]')
solvateLigand.set_parameters(molar_fractions='1.0')
solvateLigand.set_parameters(close_solvent=True)
job.add_cube(solvateLigand)

# Ligand Force Field Application
ffLigand = ForceFieldCube("ForceFieldLigand", title="Unbound Ligand Parametrization")
ffLigand.promote_parameter('ligand_forcefield', promoted_name='ligand_forcefield')
ffLigand.promote_parameter('other_forcefield', promoted_name='other_forcefield')
job.add_cube(ffLigand)

# Ligand Minimization
minimizeLigand = MDMinimizeCube("MinimizeLigand", title="Unbound Ligand Minimization")
minimizeLigand.set_parameters(restraints='noh ligand')
minimizeLigand.promote_parameter('hmr', promoted_name='hmr')
minimizeLigand.set_parameters(restraintWt=5.0)
minimizeLigand.set_parameters(center=True)
job.add_cube(minimizeLigand)

# Ligand NVT Warm-up
warmupLigand = MDNvtCube('warmupLigand', title='Unbound Ligand Warm Up')
warmupLigand.set_parameters(time=0.02)
warmupLigand.promote_parameter('temperature', promoted_name='temperature')
warmupLigand.promote_parameter('hmr', promoted_name='hmr')
warmupLigand.set_parameters(restraints="noh ligand")
warmupLigand.set_parameters(restraintWt=2.0)
warmupLigand.set_parameters(trajectory_interval=0.0)
warmupLigand.set_parameters(reporter_interval=0.0)
warmupLigand.set_parameters(suffix='warmup_ligand')
job.add_cube(warmupLigand)

# Ligand NPT Equilibration stage
equilLigand = MDNptCube('equilLigand', title='Unbound Ligand Equilibration')
equilLigand.set_parameters(time=0.02)
equilLigand.promote_parameter('temperature', promoted_name='temperature')
equilLigand.promote_parameter('pressure', promoted_name='pressure')
equilLigand.promote_parameter('hmr', promoted_name='hmr')
equilLigand.set_parameters(restraints="noh ligand")
equilLigand.set_parameters(restraintWt=0.1)
equilLigand.set_parameters(trajectory_interval=0.0)
equilLigand.set_parameters(reporter_interval=0.0)
equilLigand.set_parameters(suffix='equil')
job.add_cube(equilLigand)

sync = SyncBindingFECube("SyncCube", title="Unbound and Bound State Synchronization")
job.add_cube(sync)

report = MDFloeReportCube("report", title="Floe Report")
job.add_cube(report)

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")
job.add_cube(ofs)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")
job.add_cube(fail)


# Complex Connections
iprot.success.connect(protset.intake)
protset.success.connect(complx.protein_port)
complx.success.connect(solvateComplex.intake)
solvateComplex.success.connect(ffComplex.intake)
ffComplex.success.connect(minComplex.intake)
minComplex.success.connect(warmupComplex.intake)
warmupComplex.success.connect(equil1Complex.intake)
equil1Complex.success.connect(equil2Complex.intake)
equil2Complex.success.connect(equil3Complex.intake)
equil3Complex.success.connect(sync.intake)

# Ligand Connections
iligs.success.connect(chargelig.intake)
chargelig.success.connect(ligset.intake)
ligset.success.connect(ligid.intake)
ligid.success.connect(complx.intake)
ligid.success.connect(solvateLigand.intake)
solvateLigand.success.connect(ffLigand.intake)
ffLigand.success.connect(minimizeLigand.intake)
minimizeLigand.success.connect(warmupLigand.intake)
warmupLigand.success.connect(equilLigand.intake)
equilLigand.success.connect(sync.solvated_ligand_in_port)

# ABFE Connections
sync.success.connect(yank_proxy.intake)
yank_proxy.success.connect(report.intake)
yank_proxy.failure.connect(fail.intake)
yank_proxy.cycle_out_port.connect(abfe.intake)
report.success.connect(ofs.intake)
report.failure.connect(fail.intake)
abfe.success.connect(yank_proxy.cycle_in_port)
abfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
