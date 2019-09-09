#############################################################################
# Copyright (C) 2019 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np

from openeye import oechem

from openeye import oezap

from simtk import (unit, openmm)

from oeommtools.utils import oemol_to_openmmTop

import parmed



def ParmedAndOEMolAtomsMatch( parmedObj, mol):
    # first verify that the total number of atoms match
    pmdAtoms = parmedObj.atoms
    if len(pmdAtoms)!=mol.NumAtoms():
        print('Unequal number of atoms between parmed ({}) and OEMol({})'.
             format(len(pmdAtoms),mol.NumAtoms()))
        return False
    # second verify each atom's atomic numbers match
    for pmdAtm, oeAtom in zip(pmdAtoms,mol.GetAtoms()):
        if pmdAtm.atomic_number!=oeAtom.GetAtomicNum():
            print('different atomic number between parmed ({}) and OEMol({}) atom {}'.
                  format(pmdAtm.atomic_number,oeAtom.GetAtomicNum(),pmdAtm.name))
            return False

    return True


def ProtLigWatInteractionEFromParmedOETraj(pmed, ligOETraj, protOETraj, water_traj):

    # Generate ligand and protein parmed subsystems
    print('Generating ligand and protein parmed subsystems')

    # Parmed spitting
    pmd_split = pmed.split()

    ligand_check = False
    for idx, pmd_str in enumerate(pmd_split):

        if len(pmd_str[0].atoms) == ligOETraj.NumAtoms():
            # Checking Parmed ligand topology vs OE ligand topology
            for parm_at, oe_at in zip(pmd_str[0].atoms, ligOETraj.GetAtoms()):
                if parm_at.atomic_number != oe_at.GetAtomicNum():
                    raise ValueError("Atomic number mismatch between the Parmed ligand and "
                                     "the OpenEye ligand topologies: {} - {}".
                                     format(parm_at.atomic_number, oe_at.GetAtomicNum()))
            # Parmed ligand structure
            ligandPmed = pmd_str[0]

            # Ligand index in the parmed list splitting
            lig_idx = idx
            ligand_check = True
            break

    if not ligand_check:
        raise ValueError("The ligand cannot be found")

    atom_count = 0
    for idx in range(0, lig_idx):
        atom_count += len(pmd_split[idx][0].atoms) * len(pmd_split[idx][1])

    proteinPmed = pmed[[i for i in range(0, atom_count)]]

    # Checking Parmed Protein topology vs OE Protein topology
    if len(proteinPmed.atoms) == protOETraj.NumAtoms():

        for parm_at, oe_at in zip(proteinPmed.atoms, protOETraj.GetAtoms()):
            if parm_at.atomic_number != oe_at.GetAtomicNum():
                raise ValueError(
                    "Atomic number mismatch between the Parmed Protein and the OpenEye Protein topologies: {} - {}".
                        format(parm_at.atomic_number, oe_at.GetAtomicNum()))
    else:
        raise ValueError("Parmed Protein and OpenEye Protein number of atoms mismatch {} vs {}".
                         format(proteinPmed.atoms, protOETraj.NumAtoms()))

    # Generate Parametrized Waters
    if water_traj is not None:
        oe_water_mol = water_traj.GetActive()

        omm_wat_top, omm_wat_pos = oemol_to_openmmTop(oe_water_mol)

        wat_ff = openmm.app.ForceField("tip3p.xml")

        omm_wat_system = wat_ff.createSystem(omm_wat_top,
                                             rigidWater=False,
                                             constraints=None,
                                             nonbondedMethod=openmm.app.NoCutoff)

        water_pmed = parmed.openmm.load_topology(omm_wat_top, omm_wat_system, xyz=omm_wat_pos)

    # Parmed Complex
    complexPmed = proteinPmed + ligandPmed

    # Check that protein and ligand atoms match between parmed subsystems, ligOETraj, and protOETraj
    if not ParmedAndOEMolAtomsMatch(ligandPmed, ligOETraj):
        print('ligand atoms do not match between parmed and ligOETraj')
        return None, None, None, None
    if not ParmedAndOEMolAtomsMatch(proteinPmed, protOETraj):
        print('protein atoms do not match between parmed and protOETraj')
        return None, None, None, None
    print('protein and ligand atoms match between parmed subsystems, ligOETraj, and protOETraj')

    # While we are here put the parmed charges on protOETraj and ligOETraj as side effects
    for pmdAtom, oeAtom in zip(proteinPmed.atoms, protOETraj.GetAtoms()):
        oeAtom.SetPartialCharge(pmdAtom.charge)
    for pmdAtom, oeAtom in zip(ligandPmed.atoms, ligOETraj.GetAtoms()):
        oeAtom.SetPartialCharge(pmdAtom.charge)

    Integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin,
                                           1 / unit.picoseconds,
                                           0.002 * unit.picoseconds)

    ligOmmSys = ligandPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
    lig_integrator = openmm.LangevinIntegrator(Integrator)
    ligSim = openmm.app.Simulation(ligandPmed.topology, ligOmmSys, lig_integrator)

    # protein
    protOmmSys = proteinPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
    prot_integrator = openmm.LangevinIntegrator(Integrator)
    protSim = openmm.app.Simulation(proteinPmed.topology, protOmmSys, prot_integrator)

    # complex
    cplxOmmSys = complexPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
    cplx_integrator = openmm.LangevinIntegrator(Integrator)
    cplxSim = openmm.app.Simulation(complexPmed.topology, cplxOmmSys, cplx_integrator)

    print('OpenMM components generated for protein, ligand, and complex subsystems')

    # Water Systems
    if water_traj is not None:
        water_integrator = openmm.LangevinIntegrator(Integrator)
        water_omm_sim = openmm.app.Simulation(water_pmed.topology, omm_wat_system, water_integrator)

        lig_water_pmd = ligandPmed + water_pmed
        lig_water_sys = lig_water_pmd.createSystem(nonbondedMethod=openmm.app.NoCutoff)
        lig_water_integrator = openmm.LangevinIntegrator(Integrator)
        lig_water_sim = openmm.app.Simulation(lig_water_pmd.topology, lig_water_sys, lig_water_integrator)

        prot_water_pmd = proteinPmed + water_pmed
        prot_water_sys = prot_water_pmd.createSystem(nonbondedMethod=openmm.app.NoCutoff)
        prot_water_integrator = openmm.LangevinIntegrator(Integrator)
        prot_water_sim = openmm.app.Simulation(prot_water_pmd.topology, prot_water_sys, prot_water_integrator)

        cmplx_water_pmd = complexPmed + water_pmed
        cmplx_water_sys = cmplx_water_pmd.createSystem(nonbondedMethod=openmm.app.NoCutoff)
        cmplx_water_integrator = openmm.LangevinIntegrator(Integrator)
        cmplx_water_sim = openmm.app.Simulation(cmplx_water_pmd.topology, cmplx_water_sys, cmplx_water_integrator)

        print('OpenMM components generated for water, water-ligand, water-protein and complex-water subsystem')

    # Evaluate energies for the subsystems
    ligE = []
    protE = []
    watE = []
    cplxE = []
    plIntE = []
    lwIntE = []
    pwIntE = []
    pw_lIntE = []

    if water_traj is not None:
        zip_iter = zip(ligOETraj.GetConfs(), protOETraj.GetConfs(), water_traj.GetConfs())
    else:
        zip_iter = zip(ligOETraj.GetConfs(), protOETraj.GetConfs())

    count = 0
    for confs in zip_iter:

        ligConf = confs[0]
        protConf = confs[1]

        ligXYZdict = ligConf.GetCoords()
        ligXYZ = [openmm.Vec3(v[0], v[1], v[2]) for k, v in ligXYZdict.items()] * unit.angstrom
        ligSim.context.setPositions(ligXYZ)
        ligState = ligSim.context.getState(getEnergy=True)
        ligIntraE = ligState.getPotentialEnergy().in_units_of(
            unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
        ligE.append(ligIntraE)

        protXYZdict = protConf.GetCoords()
        protXYZ = [openmm.Vec3(v[0], v[1], v[2]) for k, v in protXYZdict.items()] * unit.angstrom
        protSim.context.setPositions(protXYZ)
        protState = protSim.context.getState(getEnergy=True)
        protIntraE = protState.getPotentialEnergy().in_units_of(
            unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
        protE.append(protIntraE)

        cplxXYZ = protXYZ + ligXYZ
        cplxSim.context.setPositions(cplxXYZ)
        cplxState = cplxSim.context.getState(getEnergy=True)
        cplxIntraE = cplxState.getPotentialEnergy().in_units_of(
            unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
        cplxE.append(cplxIntraE)

        plIntE.append(cplxIntraE - (protIntraE + ligIntraE))

        if water_traj is not None:

            water_conf = confs[2]

            # Waters intra-molecular energy
            water_xyz_dict = water_conf.GetCoords()
            water_xyz = [openmm.Vec3(v[0], v[1], v[2]) for k, v in water_xyz_dict.items()] * unit.angstrom
            water_omm_sim.context.setPositions(water_xyz)
            water_state = water_omm_sim.context.getState(getEnergy=True)
            water_intra_eng = water_state.getPotentialEnergy().in_units_of(
                unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
            watE.append(water_intra_eng)

            # Ligand-Water inter-molecular energy
            lig_water_xyz = ligXYZ + water_xyz
            lig_water_sim.context.setPositions(lig_water_xyz)
            lig_water_state = lig_water_sim.context.getState(getEnergy=True)
            lig_water_intra_eng = lig_water_state.getPotentialEnergy().in_units_of(
                unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
            lwIntE.append(lig_water_intra_eng - (ligIntraE + water_intra_eng))

            # Protein-Water inter-molecular energy
            prot_water_xyz = protXYZ + water_xyz
            prot_water_sim.context.setPositions(prot_water_xyz)
            prot_water_state = prot_water_sim.context.getState(getEnergy=True)
            prot_water_intra_eng = prot_water_state.getPotentialEnergy().in_units_of(
                unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
            pwIntE.append(prot_water_intra_eng - (protIntraE + water_intra_eng))

            cmplx_water_xyz = cplxXYZ + water_xyz
            cmplx_water_sim.context.setPositions(cmplx_water_xyz)
            cmplx_water_state = cmplx_water_sim.context.getState(getEnergy=True)
            cmplx_water_intra_eng = cmplx_water_state.getPotentialEnergy().in_units_of(
                unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole

            pw_lIntE.append(cmplx_water_intra_eng-(prot_water_intra_eng + ligIntraE))

        count += 1

    if water_traj is not None:
        print('OpenMM energies computed for protein, ligand, complex and water trajectories')
        return plIntE, cplxE, protE, ligE, watE, lwIntE, pwIntE, pw_lIntE
    else:
        print('OpenMM energies computed for protein, ligand, and complex trajectories')
        return plIntE, cplxE, protE, ligE, None, None, None, None


# def ProtLigInteractionEFromParmedOETraj( pmed, ligOETraj, protOETraj):
#     # Generate ligand and protein parmed subsystems
#     print('Generating ligand and protein parmed subsystems')
#
#     # Parmed spitting
#     pmd_split = pmed.split()
#
#     ligand_check = False
#     for idx, pmd_str in enumerate(pmd_split):
#
#         if len(pmd_str[0].atoms) == ligOETraj.NumAtoms():
#             # Checking Parmed ligand topology vs OE ligand topology
#             for parm_at, oe_at in zip(pmd_str[0].atoms, ligOETraj.GetAtoms()):
#                 if parm_at.atomic_number != oe_at.GetAtomicNum():
#                     raise ValueError("Atomic number mismatch between the Parmed ligand and "
#                                      "the OpenEye ligand topologies: {} - {}".
#                                      format(parm_at.atomic_number, oe_at.GetAtomicNum()))
#             # Parmed ligand structure
#             ligandPmed = pmd_str[0]
#
#             # Ligand index in the parmed list splitting
#             lig_idx = idx
#             ligand_check = True
#             break
#
#     if not ligand_check:
#         raise ValueError("The ligand cannot be found")
#
#     # proteinPmed = parmed.Structure()
#     #
#     # # Parmed Protein structure
#     # for idx in range(0, lig_idx):
#     #     proteinPmed += pmd_split[idx][0]
#
#     atom_count = 0
#     for idx in range(0, lig_idx):
#         atom_count += len(pmd_split[idx][0].atoms) * len(pmd_split[idx][1])
#
#     proteinPmed = pmed[[i for i in range(0, atom_count)]]
#
#     # Checking Parmed Protein topology vs OE Protein topology
#     if len(proteinPmed.atoms) == protOETraj.NumAtoms():
#
#         for parm_at, oe_at in zip(proteinPmed.atoms, protOETraj.GetAtoms()):
#             if parm_at.atomic_number != oe_at.GetAtomicNum():
#                 raise ValueError(
#                     "Atomic number mismatch between the Parmed Protein and the OpenEye Protein topologies: {} - {}".
#                         format(parm_at.atomic_number, oe_at.GetAtomicNum()))
#     else:
#         raise ValueError("Parmed Protein and OpenEye Protein number of atoms mismatch {} vs {}".
#                          format(proteinPmed.atoms, protOETraj.NumAtoms()))
#
#     # Parmed Complex
#     complexPmed = proteinPmed + ligandPmed
#
#     # Check that protein and ligand atoms match between parmed subsystems, ligOETraj, and protOETraj
#     if not ParmedAndOEMolAtomsMatch(ligandPmed, ligOETraj):
#         print('ligand atoms do not match between parmed and ligOETraj')
#         return None, None, None, None
#     if not ParmedAndOEMolAtomsMatch( proteinPmed, protOETraj):
#         print('protein atoms do not match between parmed and protOETraj')
#         return None, None, None, None
#     print('protein and ligand atoms match between parmed subsystems, ligOETraj, and protOETraj')
#     #
#     # While we are here put the parmed charges on protOETraj and ligOETraj as side effects
#     for pmdAtom, oeAtom in zip( proteinPmed.atoms, protOETraj.GetAtoms()):
#         oeAtom.SetPartialCharge( pmdAtom.charge)
#     for pmdAtom, oeAtom in zip( ligandPmed.atoms, ligOETraj.GetAtoms()):
#         oeAtom.SetPartialCharge( pmdAtom.charge)
#     #
#     # Generate openmm components needed for energy evals
#     # ligand
#     ligOmmSys = ligandPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
#     ligIntegrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
#                                               0.002 * unit.picoseconds)
#     ligSim = openmm.app.Simulation(ligandPmed.topology, ligOmmSys, ligIntegrator)
#     # protein
#     protOmmSys = proteinPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
#     protIntegrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
#                                               0.002 * unit.picoseconds)
#     protSim = openmm.app.Simulation(proteinPmed.topology, protOmmSys, protIntegrator)
#     # complex
#     cplxOmmSys = complexPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
#     cplxIntegrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
#                                               0.002 * unit.picoseconds)
#     cplxSim = openmm.app.Simulation(complexPmed.topology, cplxOmmSys, cplxIntegrator)
#     print('OpenMM components generated for protein, ligand, and complex subsystems')
#     #
#     # Evaluate energies for the subsystems
#     ligE = []
#     protE = []
#     cplxE = []
#     plIntE = []
#     #
#     for ligConf, protConf in zip(ligOETraj.GetConfs(), protOETraj.GetConfs()):
#         ligXYZdict = ligConf.GetCoords()
#         ligXYZ = [openmm.Vec3(v[0], v[1], v[2]) for k, v in ligXYZdict.items()] * unit.angstrom
#         ligSim.context.setPositions(ligXYZ)
#         ligState = ligSim.context.getState(getEnergy=True)
#         ligIntraE = ligState.getPotentialEnergy().in_units_of(
#             unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
#         ligE.append( ligIntraE)
#
#         protXYZdict = protConf.GetCoords()
#         protXYZ = [openmm.Vec3(v[0], v[1], v[2]) for k, v in protXYZdict.items()] * unit.angstrom
#         protSim.context.setPositions(protXYZ)
#         protState = protSim.context.getState(getEnergy=True)
#         protIntraE = protState.getPotentialEnergy().in_units_of(
#             unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
#         protE.append( protIntraE)
#
#         cplxXYZ = protXYZ + ligXYZ
#         cplxSim.context.setPositions(cplxXYZ)
#         cplxState = cplxSim.context.getState(getEnergy=True)
#         cplxIntraE = cplxState.getPotentialEnergy().in_units_of(
#             unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
#         cplxE.append(cplxIntraE)
#
#         plIntE.append(cplxIntraE-(protIntraE+ligIntraE))
#
#     print('OpenMM energies computed for protein, ligand, and complex trajectories')
#     return plIntE, cplxE, protE, ligE


def PBSA(ligand, protein):
    bind = oezap.OEBind()
    bind.SetProtein(protein)
    results = oezap.OEBindResults()
    if not bind.Bind(ligand, results):
        print('zap Bind run failed for {} with {}'.format(ligand.GetTitle(), protein.GetTitle()))
        return None, None, None, None, None, None
    # convert key values from kT to kcal/mol
    kTtoKcal = 0.59
    Ebind = kTtoKcal*results.GetBindingEnergy()
    EbindEl = kTtoKcal*results.GetZapEnergy()
    EdesolEl = kTtoKcal*results.GetDesolvationEnergy()
    EintEl = kTtoKcal*( results.GetZapEnergy()-results.GetDesolvationEnergy() )
    EbindSA = kTtoKcal*results.GetBuriedAreaEnergy()
    SAburied = results.GetBuriedArea()
    #
    return Ebind, EbindEl, EdesolEl, EintEl, EbindSA, SAburied


def TrajPBSA( ligOETraj, protOETraj, radiiType=oechem.OERadiiType_Zap9):
    # set ZAP9 radii on protein and ligand
    oechem.OEAssignRadii(protOETraj, radiiType, oechem.OERadiiType_HonigIonicCavity)
    oechem.OEAssignRadii(ligOETraj, radiiType)
    #
    zapBind = []
    zapBindEl = []
    zapDesolEl = []
    zapIntEl = []
    zapBindSA25 = []
    saBuried = []
    #
    for protConf, ligConf in zip(protOETraj.GetConfs(), ligOETraj.GetConfs()):
        Ebind, EbindEl, EdesolEl, EintEl, EbindSA, buriedSA = PBSA(ligConf, protConf)
        zapBind.append(Ebind)
        zapBindEl.append(EbindEl)
        zapDesolEl.append(EdesolEl)
        zapIntEl.append(EintEl)
        zapBindSA25.append(EbindSA)
        saBuried.append(buriedSA)
    #
    return zapBind, zapBindEl, zapDesolEl, zapIntEl, zapBindSA25, saBuried
