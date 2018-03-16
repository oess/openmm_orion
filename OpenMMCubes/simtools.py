import sys, mdtraj, tarfile, os
import numpy as np
from sys import stdout
from openeye import oechem
from simtk import unit, openmm
from simtk.openmm import app
from floe.api.orion import in_orion,  upload_file
from oeommtools import utils as oeommutils


def simulation(mdData, **opt):
    """
    This supporting function performs: OpenMM Minimization, NVT and NPT
    Molecular Dynamics (MD) simulations

    Parameters
    ----------
    mdData : MDData data object
        The object which recovers the relevant Parmed structure data
        to perform MD
    opt: python dictionary
        A dictionary containing all the MD setting info
    """
    
    if opt['Logger'] is None:
        printfile = sys.stdout
    else:
        printfile = opt['Logger'].file

    # MD data extracted from Parmed
    structure = mdData.structure
    topology = mdData.topology
    positions = mdData.positions
    velocities = mdData.velocities
    box = mdData.box

    # Time step in ps
    if opt['hmr']:
        stepLen = 0.004 * unit.picoseconds
        opt['Logger'].info("Hydrogen Mass reduction is On")
    else:
        stepLen = 0.002 * unit.picoseconds

    opt['timestep'] = stepLen

    # Centering the system to the OpenMM Unit Cell
    if opt['center'] and box is not None:
        opt['Logger'].info("Centering is On")
        # Numpy array in A
        coords = structure.coordinates
        # System Center of Geometry
        cog = np.mean(coords, axis=0)
        # System box vectors
        box_v = structure.box_vectors.in_units_of(unit.angstrom)/unit.angstrom
        box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])
        # Translation vector
        delta = box_v/2 - cog
        # New Coordinates
        new_coords = coords + delta
        structure.coordinates = new_coords
        positions = structure.positions

    # OpenMM system
    if box is not None:
        system = structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                        nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                        constraints=eval("app.%s" % opt['constraints']),
                                        removeCMMotion=False, hydrogenMass=4.0*unit.amu if opt['hmr'] else None)
    else:  # Vacuum
        system = structure.createSystem(nonbondedMethod=app.NoCutoff,
                                        constraints=eval("app.%s" % opt['constraints']),
                                        removeCMMotion=False, hydrogenMass=4.0*unit.amu if opt['hmr'] else None)

    # OpenMM Integrator
    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, stepLen)

    if opt['SimType'] == 'npt':
        if box is None:
            oechem.OEThrow.Fatal("NPT simulation without box vector")

        # Add Force Barostat to the system
        system.addForce(openmm.MonteCarloBarostat(opt['pressure']*unit.atmospheres, opt['temperature']*unit.kelvin, 25))

    # Apply restraints
    if opt['restraints']:
        opt['Logger'].info("RESTRAINT mask applied to: {}"
                           "\tRestraint weight: {}".format(opt['restraints'],
                                                           opt['restraintWt'] *
                                                           unit.kilocalories_per_mole/unit.angstroms**2))
        # Select atom to restraint
        res_atom_set = oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['restraints'])
        opt['Logger'].info("Number of restraint atoms: {}".format(len(res_atom_set)))
        # define the custom force to restrain atoms to their starting positions
        force_restr = openmm.CustomExternalForce('k_restr*periodicdistance(x, y, z, x0, y0, z0)^2')
        # Add the restraint weight as a global parameter in kcal/mol/A^2
        force_restr.addGlobalParameter("k_restr", opt['restraintWt']*unit.kilocalories_per_mole/unit.angstroms**2)
        # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
        force_restr.addPerParticleParameter("x0")
        force_restr.addPerParticleParameter("y0")
        force_restr.addPerParticleParameter("z0")

        for idx in range(0, len(positions)):
            if idx in res_atom_set:
                xyz = positions[idx].in_units_of(unit.nanometers)/unit.nanometers
                force_restr.addParticle(idx, xyz)
        
        system.addForce(force_restr)

    # Freeze atoms
    if opt['freeze']:
        opt['Logger'].info("FREEZE mask applied to: {}".format(opt['freeze']))

        freeze_atom_set = oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['freeze'])
        opt['Logger'].info("Number of frozen atoms: {}".format(len(freeze_atom_set)))
        # Set atom masses to zero
        for idx in range(0, len(positions)):
            if idx in freeze_atom_set:
                system.setParticleMass(idx, 0.0)

    # Platform Selection
    if opt['platform'] == 'Auto':
        # simulation = app.Simulation(topology, system, integrator)
        # Select the platform
        for plt_name in ['CUDA', 'OpenCL', 'CPU', 'Reference']:
            try:
                platform = openmm.Platform_getPlatformByName(plt_name)
                break
            except:
                if plt_name == 'Reference':
                    oechem.OEThrow.Fatal('It was not possible to select any OpenMM Platform')
                else:
                    pass
        if platform.getName() in ['CUDA', 'OpenCL']:
            for precision in ['mixed', 'single', 'double']:
                try:
                    # Set platform precision for CUDA or OpenCL
                    properties = {'Precision': precision}

                    simulation = app.Simulation(topology, system, integrator,
                                                platform=platform,
                                                platformProperties=properties)
                    break
                except:
                    if precision == 'double':
                        oechem.OEThrow.Fatal('It was not possible to select any Precision '
                                             'for the selected Platform: {}'.format(platform.getName()))
                    else:
                        pass
        else:  # CPU or Reference
            simulation = app.Simulation(topology, system, integrator, platform=platform)
    else:  # Not Auto Platform selection
        try:
            platform = openmm.Platform.getPlatformByName(opt['platform'])
        except Exception as e:
            oechem.OEThrow.Fatal('The selected platform is not supported: {}'.format(str(e)))

        if opt['platform'] in ['CUDA', 'OpenCL']:
            try:
                # Set platform CUDA or OpenCL precision
                properties = {'Precision': opt['cuda_opencl_precision']}

                simulation = app.Simulation(topology, system, integrator,
                                            platform=platform,
                                            platformProperties=properties)
            except Exception:
                oechem.OEThrow.Fatal('It was not possible to set the {} precision for the {} platform'
                                     .format(opt['cuda_opencl_precision'], opt['platform']))
        else:  # CPU or Reference Platform
            simulation = app.Simulation(topology, system, integrator, platform=platform)

    # Set starting positions and velocities
    simulation.context.setPositions(positions)

    # Set Box dimensions
    if box is not None:
        simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

    # If the velocities are not present in the Parmed structure
    # new velocity vectors are generated otherwise the system is
    # restarted from the previous State
    if opt['SimType'] in ['nvt', 'npt']:

        if opt['trajectory_interval']:
            structure.save(opt['outfname']+'.pdb', overwrite=True)
            # GAC ADDED - TESTING
            # Preserve original pdb file residue numbers
            pdbfname_test = opt['outfname'] + '_ordering_test' + '.pdb'
            ofs = oechem.oemolostream(pdbfname_test)
            flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
            ofs.SetFlavor(oechem.OEFormat_PDB, flavor)

            new_temp_mol = oeommutils.openmmTop_to_oemol(structure.topology, structure.positions, verbose=False)
            new_pos = new_temp_mol.GetCoords()
            opt['molecule'].SetCoords(new_pos)
            oechem.OEWriteConstMolecule(ofs, opt['molecule'])

        if velocities is not None:
            opt['Logger'].info('RESTARTING simulation from a previous State')
            simulation.context.setVelocities(velocities)
        else:
            # Set the velocities drawing from the Boltzmann distribution at the selected temperature
            opt['Logger'].info('GENERATING a new starting State')
            simulation.context.setVelocitiesToTemperature(opt['temperature']*unit.kelvin)

        # Convert simulation time in steps
        opt['steps'] = int(round(opt['time']/(stepLen.in_units_of(unit.picoseconds)/unit.picoseconds)))
        
        # Set Reporters
        for rep in getReporters(**opt):
            simulation.reporters.append(rep)
            
    # OpenMM platform information
    mmver = openmm.version.version
    mmplat = simulation.context.getPlatform()

    if opt['verbose']:
        # Host information
        from platform import uname
        for k, v in uname()._asdict().items():
            print(k, ':', v, file=printfile)
        # Platform properties
        for prop in mmplat.getPropertyNames():
            val = mmplat.getPropertyValue(simulation.context, prop)
            print(prop, ':', val, file=printfile)
            
    print('OpenMM({}) simulation generated for {} platform'.format(mmver, mmplat.getName()), file=printfile)

    if opt['SimType'] in ['nvt', 'npt']:

        opt['Logger'].info('Running {time} ps = {steps} steps of {SimType} at {temperature} K'.format(**opt))
        
        # Start Simulation
        simulation.step(opt['steps'])

        if box is not None:
            state = simulation.context.getState(getPositions=True, getVelocities=True,
                                                getEnergy=True, enforcePeriodicBox=True)
        else:
            state = simulation.context.getState(getPositions=True, getVelocities=True,
                                                getEnergy=True, enforcePeriodicBox=False)
        
    elif opt['SimType'] == 'min':
        
        # Run a first minimization on the Reference platform
        platform_reference = openmm.Platform.getPlatformByName('Reference')
        integrator_reference = openmm.LangevinIntegrator(opt['temperature'] * unit.kelvin,
                                                         1 / unit.picoseconds, stepLen)
        simulation_reference = app.Simulation(topology, system, integrator_reference, platform=platform_reference)
        # Set starting positions and velocities
        simulation_reference.context.setPositions(positions)

        state_reference_start = simulation_reference.context.getState(getEnergy=True)

        # Set Box dimensions
        if box is not None:
            simulation_reference.context.setPeriodicBoxVectors(box[0], box[1], box[2])

        simulation_reference.minimizeEnergy(tolerance=1e5*unit.kilojoule_per_mole)

        state_reference_end = simulation_reference.context.getState(getPositions=True)

        # Start minimization on the selected platform
        opt['Logger'].info('Minimization steps: {steps}'.format(**opt))

        # Set positions after minimization on the Reference Platform
        simulation.context.setPositions(state_reference_end.getPositions())

        simulation.minimizeEnergy(maxIterations=opt['steps'])

        state = simulation.context.getState(getPositions=True, getEnergy=True)
        print('Initial Energy = {}\nMinimized energy = {}'.format(
            state_reference_start.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole),
            state.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)),
              file=printfile)

    # OpenMM Quantity object
    structure.positions = state.getPositions(asNumpy=False)
    # OpenMM Quantity object
    if box is not None:
        structure.box_vectors = state.getPeriodicBoxVectors()

    if opt['SimType'] in ['nvt', 'npt']:
        # numpy array in units of angstrom/picoseconds
        structure.velocities = state.getVelocities(asNumpy=False)

        # If required uploading files to Orion
        _file_processing(**opt)

    # Update the OEMol complex positions to match the new
    # Parmed structure after the simulation
    new_temp_mol = oeommutils.openmmTop_to_oemol(structure.topology, structure.positions, verbose=False)
    new_pos = new_temp_mol.GetCoords()
    opt['molecule'].SetCoords(new_pos)

    return


def _file_processing(**opt):
    """
    This supporting function compresses the produced trajectory
    and supporting files in a .tar file (if required ) and eventually
    uploaded them to Orion. If not .tar file is selected then all the
    generated files are eventually uploaded in Orion

    Parameters
    ----------
    opt: python dictionary
        A dictionary containing all the MD setting info
    """

    # Set the trajectory file name
    if opt['trajectory_filetype'] == 'NetCDF':
        trj_fn = opt['outfname'] +'.nc'
    elif opt['trajectory_filetype'] == 'DCD':
        trj_fn = opt['outfname'] +'.dcd'
    elif opt['trajectory_filetype'] == 'HDF5':
        trj_fn = opt['outfname'] + '.hdf5'
    else:
        oechem.OEThrow.Fatal("The selected trajectory filetype is not supported: {}"
                             .format(opt['trajectory_filetype']))
    # Set .pdb file names
    pdb_fn = opt['outfname'] + '.pdb'
    pdb_order_fn = opt['outfname'] + '_ordering_test' + '.pdb'
    log_fn = opt['outfname'] + '.log'

    # List all the file names
    fnames = [trj_fn, pdb_fn, pdb_order_fn, log_fn]

    ex_files = []

    # Check which file names are actually produced files
    for fn in fnames:
        if os.path.isfile(fn):
            ex_files.append(fn)

    # Tar the outputted files if required
    if opt['tar']:

        tarname = opt['outfname'] + '.tar'

        opt['Logger'].info('Creating tar file: {}'.format(tarname))

        tar = tarfile.open(tarname, "w")

        for name in ex_files:
            opt['Logger'].info('Adding {} to {}'.format(name, tarname))
            tar.add(name)
        tar.close()

        opt['molecule'].SetData(oechem.OEGetTag("Tar_fname"), tarname)

        if in_orion():
            upload_file(tarname, tarname, tags=['TRJ_INFO'])

        # Clean up files that have been added to tar.
        for tmp in ex_files:
            try:
                os.remove(tmp)
            except:
                pass
    else:  # If not .tar file is required the files are eventually uploaded in Orion
        if in_orion():
            for fn in ex_files:
                upload_file(fn, fn, tags=['TRJ_INFO'])

    return


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
        (2) traj_reporter: writes trajectory to file. Supported format .nc, .dcd, .hdf5
    """
    if totalSteps is None:
        totalSteps = opt['steps']
    if outfname is None:
        outfname = opt['outfname']

    reporters = []

    if opt['reporter_interval']:

        reporter_steps = int(round(opt['reporter_interval']/(
                opt['timestep'].in_units_of(unit.picoseconds)/unit.picoseconds)))

        state_reporter = app.StateDataReporter(outfname+'.log', separator="\t",
                                               reportInterval=reporter_steps,
                                               step=True,
                                               potentialEnergy=True, totalEnergy=True,
                                               volume=True, density=True, temperature=True)

        reporters.append(state_reporter)

        progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                                  reportInterval=reporter_steps,
                                                  step=True, totalSteps=totalSteps,
                                                  time=True, speed=True, progress=True,
                                                  elapsedTime=True, remainingTime=True)

        reporters.append(progress_reporter)

    if opt['trajectory_interval']:

        trajectory_steps = int(round(opt['trajectory_interval'] / (
                opt['timestep'].in_units_of(unit.picoseconds) / unit.picoseconds)))

        trj_fname = outfname
        # Trajectory file format selection
        if opt['trajectory_filetype'] == 'NetCDF':
            trj_fname += '.nc'
            traj_reporter = mdtraj.reporters.NetCDFReporter(trj_fname, trajectory_steps)
        elif opt['trajectory_filetype'] == 'DCD':
            trj_fname += '.dcd'
            traj_reporter = app.DCDReporter(trj_fname, trajectory_steps)
        elif opt['trajectory_filetype'] == 'HDF5':
            trj_fname += '.hdf5'
            mdtraj.reporters.HDF5Reporter(trj_fname, trajectory_steps)
        else:
            oechem.OEThrow.Fatal("The selected trajectory file format is not supported: {}"
                                 .format(opt['trajectory_filetype']))

        opt['molecule'].SetData(oechem.OEGetTag("Trj_fname"), trj_fname)

        reporters.append(traj_reporter)

    return reporters