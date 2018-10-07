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


import mdtraj

import numpy as np

from sys import stdout

from simtk import (unit,
                   openmm)

from simtk.openmm import app

from oeommtools import utils as oeommutils

from platform import uname

import os

import fcntl

import time

from floe.api.orion import in_orion

from MDCubes.utils import MDState

try:
    import bz2
    have_bz2 = True
except: have_bz2 = False

try:
    import gzip
    have_gzip = True
except: have_gzip = False

import simtk.openmm as mm
import math
import time


def local_cluster(sim):
    def wrapper(*args):
        if 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(',')
            mdData = args[0]
            opt = args[1]
            opt['Logger'].warn("LOCAL FLOE CLUSTER OPTION IN USE")
            while True:

                gpu_id = gpus_available_indexes[opt['system_id'] % len(gpus_available_indexes)]

                try:
                    # fn = os.path.join('/tmp/', gpu_id + '.txt')
                    fn = gpu_id + '.txt'
                    with open(fn, 'a') as file:
                        fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        file.write("name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(
                                                                                opt['system_title'],
                                                                                opt['system_id'],
                                                                                gpus_available_indexes,
                                                                                gpu_id))
                        opt['gpu_id'] = gpu_id
                        sim(mdData, opt)
                        time.sleep(2.0)
                        fcntl.flock(file, fcntl.LOCK_UN)
                        break
                except BlockingIOError:
                    time.sleep(0.1)
                except Exception as e:  # If the simulation fails for other reasons
                    try:
                        fcntl.flock(file, fcntl.LOCK_UN)
                    except:
                        pass
                    raise ValueError("{} Simulation Failed".format(e.message))
        else:
            sim(*args)

    return wrapper


@local_cluster
def simulation(parmed_structure, opt):
    """
    This supporting function performs: OpenMM Minimization, NVT and NPT
    Molecular Dynamics (MD) simulations

    Parameters
    ----------
    parmed_structure : Parmed Structure data object
        The Object is used to recover the MD data
    opt: python dictionary
        A dictionary containing all the MD setting info
    """
    # MD data extracted from Parmed
    mdData = MDState(parmed_structure)

    topology = parmed_structure.topology
    positions = mdData.get_positions()
    velocities = mdData.get_velocities()
    box = mdData.get_box_vectors()

    # Time step in ps
    if opt['hmr']:
        stepLen = 0.004 * unit.picoseconds
        opt['Logger'].info("Hydrogen Mass reduction is On")
    else:
        stepLen = 0.002 * unit.picoseconds

    opt['timestep'] = stepLen

    # Centering the system to the OpenMM Unit Cell
    if opt['center'] and box is not None:
        opt['Logger'].info("[{}] Centering is On".format(opt['CubeTitle']))
        # Numpy array in A
        coords = parmed_structure.coordinates
        # System Center of Geometry
        cog = np.mean(coords, axis=0)
        # System box vectors
        box_v = parmed_structure.box_vectors.in_units_of(unit.angstrom)/unit.angstrom
        box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])
        # Translation vector
        delta = box_v/2 - cog
        # New Coordinates
        new_coords = coords + delta
        parmed_structure.coordinates = new_coords
        positions = parmed_structure.positions

    # OpenMM system
    if box is not None:
        system = parmed_structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                               nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                               constraints=eval("app.%s" % opt['constraints']),
                                               removeCMMotion=False, hydrogenMass=4.0*unit.amu if opt['hmr'] else None)
    else:  # Vacuum
        system = parmed_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                               constraints=eval("app.%s" % opt['constraints']),
                                               removeCMMotion=False, hydrogenMass=4.0*unit.amu if opt['hmr'] else None)

    # OpenMM Integrator
    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, stepLen)

    if opt['SimType'] == 'npt':
        if box is None:
            raise ValueError("NPT simulation without box vector")

        # Add Force Barostat to the system
        system.addForce(openmm.MonteCarloBarostat(opt['pressure']*unit.atmospheres, opt['temperature']*unit.kelvin, 25))

    # Apply restraints
    if opt['restraints']:
        opt['Logger'].info("[{}] RESTRAINT mask applied to: {}"
                           "\tRestraint weight: {}".format(opt['CubeTitle'],
                                                           opt['restraints'],
                                                           opt['restraintWt'] *
                                                           unit.kilocalories_per_mole/unit.angstroms**2))
        # Select atom to restraint
        res_atom_set = oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['restraints'])
        opt['Logger'].info("[{}] Number of restraint atoms: {}".format(opt['CubeTitle'],
                                                                       len(res_atom_set)))
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
        opt['Logger'].info("[{}] FREEZE mask applied to: {}".format(opt['CubeTitle'],
                                                                    opt['freeze']))

        freeze_atom_set = oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['freeze'])
        opt['Logger'].info("[{}] Number of frozen atoms: {}".format(opt['CubeTitle'],
                                                                    len(freeze_atom_set)))
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
                    raise ValueError('It was not possible to select any OpenMM Platform')
                else:
                    pass
        if platform.getName() in ['CUDA', 'OpenCL']:
            for precision in ['mixed', 'single', 'double']:
                try:
                    # Set platform precision for CUDA or OpenCL
                    properties = {'Precision': precision}

                    if 'gpu_id' in opt and 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():
                        properties['DeviceIndex'] = opt['gpu_id']

                    simulation = app.Simulation(topology, system, integrator,
                                                platform=platform,
                                                platformProperties=properties)
                    break
                except:
                    if precision == 'double':
                        raise ValueError('It was not possible to select any Precision '
                                         'for the selected Platform: {}'.format(platform.getName()))
                    else:
                        pass
        else:  # CPU or Reference
            simulation = app.Simulation(topology, system, integrator, platform=platform)
    else:  # Not Auto Platform selection
        try:
            platform = openmm.Platform.getPlatformByName(opt['platform'])
        except Exception as e:
            raise ValueError('The selected platform is not supported: {}'.format(str(e)))

        if opt['platform'] in ['CUDA', 'OpenCL']:
            try:
                # Set platform CUDA or OpenCL precision
                properties = {'Precision': opt['cuda_opencl_precision']}

                simulation = app.Simulation(topology, system, integrator,
                                            platform=platform,
                                            platformProperties=properties)
            except Exception:
                raise ValueError('It was not possible to set the {} precision for the {} platform'
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

        if velocities is not None:
            opt['Logger'].info('[{}] RESTARTING simulation from a previous State'.format(opt['CubeTitle']))
            simulation.context.setVelocities(velocities)
        else:
            # Set the velocities drawing from the Boltzmann distribution at the selected temperature
            opt['Logger'].info('[{}] GENERATING a new starting State'.format(opt['CubeTitle']))
            simulation.context.setVelocitiesToTemperature(opt['temperature']*unit.kelvin)

        # Convert simulation time in steps
        opt['steps'] = int(round(opt['time']/(stepLen.in_units_of(unit.nanoseconds)/unit.nanoseconds)))
        
        # Set Reporters
        for rep in getReporters(**opt):
            simulation.reporters.append(rep)
            
    # OpenMM platform information
    mmver = openmm.version.version
    mmplat = simulation.context.getPlatform()

    str_logger = '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32
    str_logger += '\n' + '{:<25} = {:<10}'.format('time step', str(opt['timestep']))

    # Host information
    for k, v in uname()._asdict().items():
        str_logger += "\n{:<25} = {:<10}".format(k, v)
        opt['Logger'].info("[{}] {} : {}".format(opt['CubeTitle'],
                                                 k, v))

    # Platform properties
    for prop in mmplat.getPropertyNames():
        val = mmplat.getPropertyValue(simulation.context, prop)
        str_logger += "\n{:<25} = {:<10}".format(prop, val)
        opt['Logger'].info("[{}] {} : {}".format(opt['CubeTitle'],
                                                 prop, val))

    info = "{:<25} = {:<10}".format("OpenMM Version", mmver)
    opt['Logger'].info("[{}] OpenMM Version : {}".format(opt['CubeTitle'], mmver))
    str_logger += '\n'+info

    info = "{:<25} = {:<10}".format("Platform in use", mmplat.getName())
    opt['Logger'].info("[{}] Platform in use : {}".format(opt['CubeTitle'], mmplat.getName()))
    str_logger += '\n'+info

    if opt['SimType'] in ['nvt', 'npt']:

        if opt['SimType'] == 'nvt':
            opt['Logger'].info("[{}] Running time : {time} ns => {steps} steps of {SimType} at "
                               "{temperature} K".format(opt['CubeTitle'], **opt))
            info = "{:<25} = {time} ns => {steps} steps of {SimType} at " \
                   "{temperature} K".format("Running time", **opt)
        else:
            opt['Logger'].info(
                "[{}] Running time : {time} ns => {steps} steps of {SimType} "
                "at {temperature} K pressure {pressure} atm".format(opt['CubeTitle'], **opt))
            info = "{:<25} = {time} ns => {steps} steps of {SimType} at " \
                   "{temperature} K pressure {pressure} atm".format("Running time", **opt)

        str_logger += '\n' + info

        if opt['trajectory_interval']:

            total_frames = int(round(opt['time'] / opt['trajectory_interval']))

            opt['Logger'].info('[{}] Total trajectory frames : {}'.format(opt['CubeTitle'], total_frames))
            info = '{:<25} = {:<10}'.format('Total trajectory frames', total_frames)
            str_logger += '\n' + info

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
        if opt['steps'] == 0:
            opt['Logger'].info('[{}] Minimization steps: until convergence is found'.format(opt['CubeTitle']))
        else:
            opt['Logger'].info('[{}] Minimization steps: {steps}'.format(opt['CubeTitle'], **opt))

        # Set positions after minimization on the Reference Platform
        simulation.context.setPositions(state_reference_end.getPositions())

        simulation.minimizeEnergy(maxIterations=opt['steps'])

        state = simulation.context.getState(getPositions=True, getEnergy=True)

        ie = '{:<25} = {:<10}'.format('Initial Potential Energy', str(state_reference_start.getPotentialEnergy().
                                                                      in_units_of(unit.kilocalorie_per_mole)))
        fe = '{:<25} = {:<10}'.format('Minimized Potential Energy', str(state.getPotentialEnergy().
                                                                        in_units_of(unit.kilocalorie_per_mole)))

        opt['Logger'].info(ie)
        opt['Logger'].info(fe)

        str_logger += '\n' + ie + '\n' + fe

    # OpenMM Quantity object
    parmed_structure.positions = state.getPositions(asNumpy=False)
    # OpenMM Quantity object
    if box is not None:
        parmed_structure.box_vectors = state.getPeriodicBoxVectors()

    if opt['SimType'] in ['nvt', 'npt']:
        # numpy array in units of angstrom/picoseconds
        parmed_structure.velocities = state.getVelocities(asNumpy=False)

    # Update the OEMol complex positions to match the new
    # Parmed structure after the simulation
    new_temp_mol = oeommutils.openmmTop_to_oemol(parmed_structure.topology, parmed_structure.positions, verbose=False)
    new_pos = new_temp_mol.GetCoords()
    opt['molecule'].SetCoords(new_pos)

    # Update the string logger
    opt['str_logger'] += str_logger

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
                opt['timestep'].in_units_of(unit.nanoseconds)/unit.nanoseconds)))

        state_reporter = app.StateDataReporter(outfname+'.log', separator="\t",
                                               reportInterval=reporter_steps,
                                               step=True,
                                               potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                               volume=True, density=True, temperature=True)

        reporters.append(state_reporter)

        progress_reporter = StateDataReporterName(stdout, system_name=opt['system_title']+'_'+str(opt['system_id']),
                                                  separator="\t",
                                                  reportInterval=reporter_steps,
                                                  step=False, totalSteps=totalSteps,
                                                  time=True, speed=True, progress=True,
                                                  elapsedTime=False, remainingTime=True)

        reporters.append(progress_reporter)

    if opt['trajectory_interval']:

        trajectory_steps = int(round(opt['trajectory_interval'] / (
                opt['timestep'].in_units_of(unit.nanoseconds) / unit.nanoseconds)))

        traj_reporter = mdtraj.reporters.HDF5Reporter(outfname+'.h5', trajectory_steps)

        reporters.append(traj_reporter)

    return reporters


class StateDataReporterName(object):
    """
    This class has been adapted From OpenMM 7.1.1 to print the system name that is in process

    StateDataReporter outputs information about a simulation, such as energy and temperature, to a file.

    To use it, create a StateDataReporter, then add it to the Simulation's list of reporters.  The set of
    data to write is configurable using boolean flags passed to the constructor.  By default the data is
    written in comma-separated-value (CSV) format, but you can specify a different separator to use.
    """

    def __init__(self, file, reportInterval, system_name=None, step=False,
                 time=False, potentialEnergy=False, kineticEnergy=False,
                 totalEnergy=False, temperature=False, volume=False, density=False,
                 progress=False, remainingTime=False, speed=False, elapsedTime=False,
                 separator=',', systemMass=None, totalSteps=None):
        """Create a StateDataReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reportInterval : int
            The interval (in time steps) at which to write frames
        system_name : string=None
            The string that specify the system name
        step : bool=False
            Whether to write the current step index to the file
        time : bool=False
            Whether to write the current time to the file
        potentialEnergy : bool=False
            Whether to write the potential energy to the file
        kineticEnergy : bool=False
            Whether to write the kinetic energy to the file
        totalEnergy : bool=False
            Whether to write the total energy to the file
        temperature : bool=False
            Whether to write the instantaneous temperature to the file
        volume : bool=False
            Whether to write the periodic box volume to the file
        density : bool=False
            Whether to write the system density to the file
        progress : bool=False
            Whether to write current progress (percent completion) to the file.
            If this is True, you must also specify totalSteps.
        remainingTime : bool=False
            Whether to write an estimate of the remaining clock time until
            completion to the file.  If this is True, you must also specify
            totalSteps.
        speed : bool=False
            Whether to write an estimate of the simulation speed in ns/day to
            the file
        elapsedTime : bool=False
            Whether to write the elapsed time of the simulation in seconds to
            the file.
        separator : string=','
            The separator to use between columns in the file
        systemMass : mass=None
            The total mass to use for the system when reporting density.  If
            this is None (the default), the system mass is computed by summing
            the masses of all particles.  This parameter is useful when the
            particle masses do not reflect their actual physical mass, such as
            when some particles have had their masses set to 0 to immobilize
            them.
        totalSteps : int=None
            The total number of steps that will be included in the simulation.
            This is required if either progress or remainingTime is set to True,
            and defines how many steps will indicate 100% completion.
        """
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if (progress or remainingTime) and totalSteps is None:
            raise ValueError('Reporting progress or remaining time requires total steps to be specified')
        if self._openedFile:
            # Detect the desired compression scheme from the filename extension
            # and open all files unbuffered
            if file.endswith('.gz'):
                if not have_gzip:
                    raise RuntimeError("Cannot write .gz file because Python could not import gzip library")
                self._out = gzip.GzipFile(fileobj=open(file, 'wb', 0))
            elif file.endswith('.bz2'):
                if not have_bz2:
                    raise RuntimeError("Cannot write .bz2 file because Python could not import bz2 library")
                self._out = bz2.BZ2File(file, 'w', 0)
            else:
                self._out = open(file, 'w')
        else:
            self._out = file
        self._system_name = system_name
        self._step = step
        self._time = time
        self._potentialEnergy = potentialEnergy
        self._kineticEnergy = kineticEnergy
        self._totalEnergy = totalEnergy
        self._temperature = temperature
        self._volume = volume
        self._density = density
        self._progress = progress
        self._remainingTime = remainingTime
        self._speed = speed
        self._elapsedTime = elapsedTime
        self._separator = separator
        self._totalMass = systemMass
        self._totalSteps = totalSteps
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = False
        self._needEnergy = potentialEnergy or kineticEnergy or totalEnergy or temperature

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, self._needsPositions, self._needsVelocities, self._needsForces, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            print('#"%s"' % ('"'+self._separator+'"').join(headers), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._initialClockTime = time.time()
            self._initialSimulationTime = state.getTime()
            self._initialSteps = simulation.currentStep
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        # Write the values.
        print(self._separator.join(str(v) for v in values), file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass

    def _constructReportValues(self, simulation, state):
        """Query the simulation for the current state of our observables of interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        Returns
        -------
        A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []
        box = state.getPeriodicBoxVectors()
        volume = box[0][0]*box[1][1]*box[2][2]
        clockTime = time.time()
        if self._system_name:
            values.append('%-10s' % self._system_name)
        if self._progress:
            values.append('%.1f%%' % (100.0*simulation.currentStep/self._totalSteps))
        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            values.append(state.getTime().value_in_unit(unit.picosecond))
        if self._potentialEnergy:
            values.append(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        if self._kineticEnergy:
            values.append(state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole))
        if self._totalEnergy:
            values.append((state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole))
        if self._temperature:
            values.append((2*state.getKineticEnergy()/(self._dof*unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin))
        if self._volume:
            values.append(volume.value_in_unit(unit.nanometer**3))
        if self._density:
            values.append((self._totalMass/volume).value_in_unit(unit.gram/unit.item/unit.milliliter))
        if self._speed:
            elapsedDays = (clockTime-self._initialClockTime)/86400.0
            elapsedNs = (state.getTime()-self._initialSimulationTime).value_in_unit(unit.nanosecond)
            if elapsedDays > 0.0:
                values.append('%.3g' % (elapsedNs/elapsedDays))
            else:
                values.append('--')
        if self._elapsedTime:
            values.append(time.time() - self._initialClockTime)
        if self._remainingTime:
            elapsedSeconds = clockTime-self._initialClockTime
            elapsedSteps = simulation.currentStep-self._initialSteps
            if elapsedSteps == 0:
                value = '--'
            else:
                estimatedTotalSeconds = (self._totalSteps-self._initialSteps)*elapsedSeconds/elapsedSteps
                remainingSeconds = int(estimatedTotalSeconds-elapsedSeconds)
                remainingDays = remainingSeconds//86400
                remainingSeconds -= remainingDays*86400
                remainingHours = remainingSeconds//3600
                remainingSeconds -= remainingHours*3600
                remainingMinutes = remainingSeconds//60
                remainingSeconds -= remainingMinutes*60
                if remainingDays > 0:
                    value = "%d:%d:%02d:%02d" % (remainingDays, remainingHours, remainingMinutes, remainingSeconds)
                elif remainingHours > 0:
                    value = "%d:%02d:%02d" % (remainingHours, remainingMinutes, remainingSeconds)
                elif remainingMinutes > 0:
                    value = "%d:%02d" % (remainingMinutes, remainingSeconds)
                else:
                    value = "0:%02d" % remainingSeconds
            values.append(value)
        return values

    def _initializeConstants(self, simulation):
        """Initialize a set of constants required for the reports

        Parameters
        - simulation (Simulation) The simulation to generate a report for
        """
        system = simulation.system
        if self._temperature:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0*unit.dalton:
                    dof += 3
            dof -= system.getNumConstraints()
            if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                dof -= 3
            self._dof = dof
        if self._density:
            if self._totalMass is None:
                # Compute the total system mass.
                self._totalMass = 0*unit.dalton
                for i in range(system.getNumParticles()):
                    self._totalMass += system.getParticleMass(i)
            elif not unit.is_quantity(self._totalMass):
                self._totalMass = self._totalMass*unit.dalton

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being reported on.
        """
        headers = []
        if self._system_name:
            headers.append('System')
        if self._progress:
            headers.append('Progress (%)')
        if self._step:
            headers.append('Step')
        if self._time:
            headers.append('Time (ps)')
        if self._potentialEnergy:
            headers.append('Potential Energy (kJ/mole)')
        if self._kineticEnergy:
            headers.append('Kinetic Energy (kJ/mole)')
        if self._totalEnergy:
            headers.append('Total Energy (kJ/mole)')
        if self._temperature:
            headers.append('Temperature (K)')
        if self._volume:
            headers.append('Box Volume (nm^3)')
        if self._density:
            headers.append('Density (g/mL)')
        if self._speed:
            headers.append('Speed (ns/day)')
        if self._elapsedTime:
            headers.append('Elapsed Time (s)')
        if self._remainingTime:
            headers.append('Time Remaining')
        return headers

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._needEnergy:
            energy = (state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole)
            if math.isnan(energy):
                raise ValueError('Energy is NaN')
            if math.isinf(energy):
                raise ValueError('Energy is infinite')

    def __del__(self):
        if self._openedFile:
            self._out.close()
