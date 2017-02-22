import io, os, base64, parmed, mdtraj, pdbfixer
import numpy as np
from sys import stdout
from openeye import oechem
from simtk import unit, openmm
from simtk.openmm import app
import OpenMMCubes.utils as utils
try:
    import cPickle as pickle
except ImportError:
    import pickle

def genProteinStructure(proteinpdb, **opt):
    #Generate protein Structure object
    forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
    protein_system = forcefield.createSystem( proteinpdb.topology )
    protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                    protein_system,
                                                    xyz=proteinpdb.positions)
    return protein_structure

def combinePostions(proteinPositions, molPositions):
    positions_unit = unit.angstroms
    positions0_dimensionless = np.array(proteinPositions / positions_unit)
    positions1_dimensionless = np.array(molPositions / positions_unit)
    coordinates = np.vstack(
        (positions0_dimensionless, positions1_dimensionless))
    natoms = len(coordinates)
    positions = np.zeros([natoms, 3], np.float32)
    for index in range(natoms):
            (x, y, z) = coordinates[index]
            positions[index, 0] = x
            positions[index, 1] = y
            positions[index, 2] = z
    positions = unit.Quantity(positions, positions_unit)
    return positions

def mergeStructure(proteinStructure, molStructure):
    structure = proteinStructure + molStructure
    positions = combinePostions(proteinStructure.positions, molStructure.positions)
    # Concatenate positions arrays (ensures same units)
    structure.positions = positions
    # Restore original box vectors
    structure.box = proteinStructure.box
    return structure

def solvateComplexStructure(structure, **opt):
    log = opt['logger']
    tmpfile = opt['outfname']+'-pl.tmp'
    structure.save(tmpfile,format='pdb',overwrite=True)

    seqres = False
    with open(tmpfile, 'r') as infile:
        for line in infile:
            if 'SEQRES' in line:
                seqres = True
                break
    if not seqres:
        log.warn('Did not find SEQRES in PDB. PDBFixer will not find missing Residues.')

    # Solvate with PDBFixer
    log.info('PDBFixer solvating {}'.format(opt['outfname']))
    log.info('\tpH = {}'.format(opt['pH']))
    log.info('\tpadding = {}'.format(unit.Quantity(opt['solvent_padding'], unit.angstroms)))
    log.info('\tionicStrength = {}'.format(unit.Quantity(opt['salt_concentration'], unit.millimolar)))
    fixer = pdbfixer.PDBFixer(tmpfile)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()

    if fixer.missingAtoms:
        log.info('Found missing Atoms:', fixer.missingAtoms)
    if fixer.missingTerminals:
        log.info('Found missing Terminals:', fixer.missingTerminals)
    if fixer.nonstandardResidues:
        log.info('Found nonstandard Residues:', fixer.nonstandardResidues)

    fixer.replaceNonstandardResidues()
    #fixer.removeHeterogens(False)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(opt['pH'])
    fixer.addSolvent(padding=unit.Quantity(opt['solvent_padding'], unit.angstroms),
                ionicStrength=unit.Quantity(opt['salt_concentration'], unit.millimolar))

    # Load PDBFixer object back to Structure
    tmp = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)

    # Remove ligand from protein Structure by AmberMask selection
    tmp.strip(":LIG")
    tmp.save(opt['outfname']+'-nomol.tmp',format='pdb',overwrite=True)
    # Reload PDBFile
    nomol = app.PDBFile(opt['outfname']+'-nomol.tmp')
    forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
    nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)
    # Regenerate parameterized solvated protein structure
    solv_structure = parmed.openmm.load_topology(nomol.topology,
                                                nomol_system,
                                                xyz=nomol.positions)
    # Restore box vectors
    solv_structure.box = tmp.box

    tmpfiles = [ opt['outfname']+'-pl.tmp', opt['outfname']+'-nomol.tmp' ]
    utils.cleanup(tmpfiles)

    return solv_structure

def genSimFromStruct(structure, **opt):
    """ ParmedEd createSystem Defaults
    nonbondedMethod=None,
    nonbondedCutoff=Quantity(value=8.0, unit=angstrom),
    switchDistance=Quantity(value=0.0, unit=angstrom),
    constraints=None,
    rigidWater=True,
    implicitSolvent=None,
    implicitSolventKappa=None,
    implicitSolventSaltConc=Quantity(value=0.0, unit=mole/liter),
    temperature=Quantity(value=298.15, unit=kelvin),
    soluteDielectric=1.0,
    solventDielectric=78.5,
    useSASA=False,
    removeCMMotion=True,
    hydrogenMass=None,
    ewaldErrorTolerance=0.0005,
    flexibleConstraints=True,
    verbose=False,
    splitDihedrals=False
    """
    system = structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                    nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                    constraints=eval("app.%s" % opt['constraints']),
                                    temperature=opt['temperature'])

    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
    simulation = app.Simulation(structure.topology, system, integrator)
    # Set initial positions/velocities
    # Will get overwritten from saved State.
    simulation.context.setPositions(structure.positions)
    simulation.context.setVelocitiesToTemperature(opt['temperature']*unit.kelvin)
    return simulation

def minimizeSimulation(simulation, **opt):
    log = opt['logger']
    init = simulation.context.getState(getEnergy=True)
    log.info('Initial energy is {}'.format(init.getPotentialEnergy()))
    simulation.minimizeEnergy()
    st = simulation.context.getState(getPositions=True,getEnergy=True)
    log.info('Minimized energy is {}'.format(st.getPotentialEnergy()))
    return simulation

def getReporters(**opt):
    progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                        reportInterval=opt['reporter_interval'],
                                        totalSteps=opt['steps'],
                                        time=True, speed=True, progress=True,
                                        elapsedTime=True, remainingTime=True)

    state_reporter = app.StateDataReporter(opt['outfname']+'.log', separator="\t",
                                        reportInterval=opt['reporter_interval'],
                                        step=True,
                                        potentialEnergy=True, totalEnergy=True,
                                        volume=True, temperature=True)

    traj_reporter = mdtraj.reporters.HDF5Reporter(opt['outfname']+'.h5', opt['reporter_interval'])

    reporters = [state_reporter, progress_reporter, traj_reporter]
    return reporters
