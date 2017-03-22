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
    """
    Starting from OpenMM PDBFile, generates the OpenMM System for the protein and
    then create the parametrized ParmEd Structure of the protein.

    Parameters
    ----------
    proteinpdb : openmm.app.PDBFile object,
        Loaded PDBFile object of the protein.
    protein_forcefield : (opt), xml file, default='amber99sbildn.xml'
        Forcefield parameters for protein
    solvent_forcefield : opt), xml file, default='tip3p.xml'
        Forcefield parameters for solvent

    Returns
    -------
    solv_structure : parmed.structure.Structure
        The parameterized Structure of the protein with solvent molecules. (No ligand).
    """
    #Generate protein Structure object
    forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
    protein_system = forcefield.createSystem( proteinpdb.topology )
    protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                    protein_system,
                                                    xyz=proteinpdb.positions)
    return protein_structure

def combinePostions(proteinPositions, molPositions):
    """
    Loops through the positions from the ParmEd structures of the protein and ligand,
    divides by unit.angstroms which will ensure both positions arrays are in the same units.

    Parameters
    ----------
    proteinPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the protein Structure.
    molPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the molecule Structure.

    Returns
    -------
    positions : list of 3-element Quantity tuples.
        ex. unit.Quantity(positions, positions_unit)
        Combined positions of the protein and molecule Structures.
    """
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
    """
    Combines the parametrized ParmEd structures of the protein and ligand to
    create the Structure for the protein:ligand complex, while retaining the SMIRFF
    parameters on the ligand. Preserves positions and box vectors.
    (Not as easily achieved using native OpenMM tools).

    Parameters
    ----------
    proteinStructure : parmed.structure.Structure
        The parametrized structure of the protein.
    moleculeStructure : parmed.structure.Structure
        The parametrized structure of the ligand.

    Returns
    -------
    structure : parmed.structure.Structure
        The parametrized structure of the protein:ligand complex.
    """
    structure = proteinStructure + molStructure
    positions = combinePostions(proteinStructure.positions, molStructure.positions)
    # Concatenate positions arrays (ensures same units)
    structure.positions = positions
    # Restore original box vectors
    structure.box = proteinStructure.box
    return structure

def solvateComplexStructure(structure, **opt):
    """
    Uses PDBFixer to add missing atoms, assign protonation states, and solvate
    the parametrized ParmEd Structure of the protein:ligand complex. Returns the
    parametrized Structure without the ligand present.
    We initially combine the ligand with the protein to prevent PDBFixer
    from adding solvent into the binding pocket.

    Parameters
    ----------
    structure : parmed.structure.Structure
        The parametrized structure of the protein:ligand complex
    pH : (opt), float, default=7.4
        Solvent pH used to select appropriate protein protonation state.
    solvent_padding : (opt), float, default=10
        Padding around protein for solvent box (unit.angstroms)
    salt_concentration : (opt), float, default=10
        Salt concentration (millimolar)
    protein_forcefield : (opt), xml file, default='amber99sbildn9sbildn.xml'
        Forcefield parameters for protein
    solvent_forcefield : opt), xml file, default='tip3p.xml'
        Forcefield parameters for solvent

    Returns
    -------
    solv_structure : parmed.structure.Structure
        The parameterized Structure of the protein with solvent molecules. (No ligand).
    """
    tmpfile = opt['outfname']+'-pl.tmp'
    structure.save(tmpfile,format='pdb')

    seqres = False
    with open(tmpfile, 'r') as infile:
        for line in infile:
            if 'SEQRES' in line:
                seqres = True
                break
    if not seqres:
        print('Did not find SEQRES in PDB. PDBFixer will not find missing Residues.')

    # Solvate with PDBFixer
    print('PDBFixer settings for {}'.format(opt['outfname']))
    print('\tpH = {}'.format(opt['pH']))
    print('\tpadding = {}'.format(unit.Quantity(opt['solvent_padding'], unit.angstroms)))
    print('\tsalt conc. = {}'.format(unit.Quantity(opt['salt_concentration'], unit.millimolar)))
    fixer = pdbfixer.PDBFixer(tmpfile)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()

    if fixer.missingAtoms:
        print('Found missing Atoms:', fixer.missingAtoms)
    if fixer.missingTerminals:
        print('Found missing Terminals:', fixer.missingTerminals)
    if fixer.nonstandardResidues:
        print('Found nonstandard Residues:', fixer.nonstandardResidues)

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
    tmp.save(opt['outfname']+'-nomol.tmp',format='pdb')
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

def genSimFromStruct(structure, platform=None, **opt):
    """
    Uses ParmEd Structure to generate the OpenMM System,
    required to create the OpenMM Simulation. ParmEd Structures are used to keep
    the SMIRFF parameters on the molecule, but aren't necessary after storing the State.
    (Keep it anyways.)

    Parameters
    ----------
    structure : parmed.structure.Structure
        The parametrized structure of the solvated protein:ligand complex
    nonbondedMethod : (opt), str
        Can be: NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.
    nonbondedCutoff : (opt), float
        The nonbonded cutoff (in unit.angstroms), ignored if nonbondedMethod is NoCutoff.
    constraints : (opt), str
        Can be: None, HBonds, HAngles, or AllBonds.

    Returns
    -------
    simulation : openmm.app.simulation.Simulation
        The OpenMM Simulation object generated from the ParmEd Structure.
    """
    system = structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                    nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                    constraints=eval("app.%s" % opt['constraints']))
    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
    if platform is None:
        #Use the fastest available platform
        simulation = app.Simulation(structure.topology, system, integrator)
    else:
        platform=openmm.Platform.getPlatformByName(platform)
        simulation = app.Simulation(structure.topology, system, integrator, platform)

    mmver = openmm.version.version
    platform = simulation.context.getPlatform().getName()
    print('OpenMM({}) simulation generated for {} platform'.format(mmver, platform))
    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
    simulation = app.Simulation(structure.topology, system, integrator)
    # Set initial positions/velocities
    # Will get overwritten from saved State.
    simulation.context.setPositions(structure.positions)
    simulation.context.setVelocitiesToTemperature(opt['temperature']*unit.kelvin)
    return simulation

def minimizeSimulation(simulation):
    """
    Minimizes the OpenMM Simulations.

    Parameters
    ----------
    simulation : openmm.app.simulation.Simulation
        The OpenMM Simulation object generated from the ParmEd Structure to be minimized.

    Returns
    -------
    simulation : openmm.app.simulation.Simulation
        The OpenMM Simulation after minimization.
    """
    init = simulation.context.getState(getEnergy=True)
    print('Initial energy = {}'.format(init.getPotentialEnergy()))
    simulation.minimizeEnergy()
    minene = simulation.context.getState(getPositions=True,getEnergy=True).getPotentialEnergy()
    print('Minimized energy = {}'.format(minene))
    return minene, simulation

def getReporters(totalSteps=None, outfname=None, **opt):
    """
    Creates 3 OpenMM Reporters for the simulation.

    Parameters
    ----------
    totalSteps : int
        The total number of simulation steps
    reportInterval : (opt), int, default=1000
        Step frequency to write to reporter file.
    outfname : str
        Specifies the filename prefix for the reporters.

    Returns
    -------
    reporters : list of three openmm.app.simulation.reporters
        (0) state_reporter: writes energies to '.log' file.
        (1) progress_reporter: prints simulation progress to 'sys.stdout'
        (2) traj_reporter: writes trajectory to '.nc' file. AMBER NetCDF(3.0)
    """
    if totalSteps is None:
        totalSteps = opt['steps']
    if outfname is None:
        outfname = opt['outfname']
    progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                        reportInterval=opt['reporter_interval'],
                                        step=True, totalSteps=totalSteps,
                                        time=True, speed=True, progress=True,
                                        elapsedTime=True, remainingTime=True)

    state_reporter = app.StateDataReporter(outfname+'.log', separator="\t",
                                        reportInterval=opt['reporter_interval'],
                                        step=True,
                                        potentialEnergy=True, totalEnergy=True,
                                        volume=True, temperature=True)

    # Uses NetCDF(4.0) but not VMD compatible.
    #traj_reporter = mdtraj.reporters.HDF5Reporter(outfname+'.h5', opt['reporter_interval'])

    # Default to NetCDF since VMD compatible.
    traj_reporter = mdtraj.reporters.NetCDFReporter(outfname+'.nc', opt['reporter_interval'])

    reporters = [state_reporter, progress_reporter, traj_reporter]
    return reporters

def mdTrajConvert(simulation=None, outfname=None, trajectory_selection=None, trajectory_filetype='NetCDF'):
    """
    Used to convert the (mdTraj) trajectory file between: HDF5, DCD, or NetCDF.
    Can also be used to write out a subset of the trajectory.
    See mdTraj atom selection docs:
        http://mdtraj.org/1.8.0/examples/atom-selection.html

    Parameters
    ----------
    simulation : openmm.app.simulation.Simulation
        The OpenMM Simulation object used to generate the mdTraj topology for
        writing to DCD or NetCDF files.
    trajectory_selection : str, default=None (writes all atoms)
        Examples: 'protein', 'resname LIG', 'protein or resname LIG'
        See mdTraj atom selection docs:
            http://mdtraj.org/1.8.0/examples/atom-selection.html
    trajectory_filetype : str, default='NetCDF'
        Can be HDF5, DCD, or NetCDF.
    outfname : str
        Specifies the filename prefix for the trajectory files.

    """
    atom_indices = None
    traj_dict = { 'HDF5' : outfname +'.h5',
                  'DCD' : outfname +'.dcd',
                  'NetCDF' : outfname +'.nc' }
    for rep in simulation.reporters:
        try:
            trajfile = rep._traj_file
            trajfname = trajfile._handle.filename
            if trajectory_selection is not None:
                atom_indices = trajfile.topology.select('%s' % trajectory_selection)
            trajfile.close()
        except Exception as e:
            pass
    top = mdtraj.Topology.from_openmm(simulation.topology)
    outfile = traj_dict.get(trajectory_filetype)
    traj = mdtraj.load(trajfname, atom_indices=atom_indices, top=top)
    print("\tTrajectory subset: '{}'\n\t{}".format(trajectory_selection, traj))
    print("\tConverted trajectory to %s" % (outfile) )
    traj.save(outfile)
