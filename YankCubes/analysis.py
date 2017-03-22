#!/usr/local/bin/env python

# =============================================================================================
# Analyze datafiles produced by YANK.
# =============================================================================================

# =============================================================================================
# REQUIREMENTS
#
# The netcdf4-python module is now used to provide netCDF v4 support:
# http://code.google.com/p/netcdf4-python/
#
# This requires NetCDF with version 4 and multithreading support, as well as HDF5.
# =============================================================================================

import os
import os.path

import yaml
import numpy as np

import netCDF4 as netcdf  # netcdf4-python

from pymbar import MBAR  # multistate Bennett acceptance ratio
from pymbar import timeseries  # for statistical inefficiency analysis

import mdtraj
import simtk.unit as units

from . import utils

import logging
logger = logging.getLogger(__name__)

# =============================================================================================
# PARAMETERS
# =============================================================================================

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

# =============================================================================================
# SUBROUTINES
# =============================================================================================


def generate_mixing_statistics(ncfile, nequil=0):
    """
    Generate the mixing statistics

    Parameters
    ----------
    ncfile : netCDF4.Dataset
       NetCDF file
    nequil : int, optional, default=0
       If specified, only samples nequil:end will be used in analysis

    Returns
    -------
    mixing_stats : np.array of shape [nstates, nstates]
        Transition matrix estimate
    mu : np.array
        Eigenvalues of the Transition matrix sorted in descending order
    """

    # Get dimensions.
    niterations = ncfile.variables['states'].shape[0]
    nstates = ncfile.variables['states'].shape[1]

    # Compute empirical transition count matrix.
    Nij = np.zeros([nstates, nstates], np.float64)
    for iteration in range(nequil, niterations-1):
        for ireplica in range(nstates):
            istate = ncfile.variables['states'][iteration, ireplica]
            jstate = ncfile.variables['states'][iteration+1, ireplica]
            Nij[istate, jstate] += 1

    # Compute transition matrix estimate.
    # TODO: Replace with maximum likelihood reversible count estimator from msmbuilder or pyemma.
    Tij = np.zeros([nstates,nstates], np.float64)
    for istate in range(nstates):
        denom = (Nij[istate,:].sum() + Nij[:,istate].sum())
        if denom > 0:
            for jstate in range(nstates):
                Tij[istate, jstate] = (Nij[istate, jstate] + Nij[jstate, istate]) / denom
        else:
            Tij[istate, istate] = 1.0

    # Estimate eigenvalues
    mu = np.linalg.eigvals(Tij)
    mu = -np.sort(-mu)  # Sort in descending order

    return Tij, mu


def show_mixing_statistics(ncfile, cutoff=0.05, nequil=0):
    """
    Print summary of mixing statistics. Passes information off to generate_mixing_statistics then prints it out to
    the logger

    Parameters
    ----------

    ncfile : netCDF4.Dataset
       NetCDF file
    cutoff : float, optional, default=0.05
       Only transition probabilities above 'cutoff' will be printed
    nequil : int, optional, default=0
       If specified, only samples nequil:end will be used in analysis

    """

    Tij, mu = generate_mixing_statistics(ncfile, nequil=nequil)

    # Print observed transition probabilities.
    nstates = ncfile.variables['states'].shape[1]
    logger.info("Cumulative symmetrized state mixing transition matrix:")
    str_row = "%6s" % ""
    for jstate in range(nstates):
        str_row += "%6d" % jstate
    logger.info(str_row)

    for istate in range(nstates):
        str_row = ""
        str_row += "%-6d" % istate
        for jstate in range(nstates):
            P = Tij[istate,jstate]
            if P >= cutoff:
                str_row += "%6.3f" % P
            else:
                str_row += "%6s" % ""
        logger.info(str_row)

    # Estimate second eigenvalue and equilibration time.
    if mu[1] >= 1:
        logger.info("Perron eigenvalue is unity; Markov chain is decomposable.")
    else:
        logger.info("Perron eigenvalue is {0:9.5f}; state equilibration timescale is ~ {1:.1f} iterations".format(
            mu[1], 1.0 / (1.0 - mu[1]))
        )

    return


def extract_ncfile_energies(ncfile, ndiscard=0, nuse=None, g=None):
    """
    Extract and decorelate energies from the ncfile to gather common data for other functions

    Parameters
    ----------
    ncfile : NetCDF
       Input YANK netcdf file
    ndiscard : int, optional, default=0
       Number of iterations to discard to equilibration
    nuse : int, optional, default=None
       Maximum number of iterations to use (after discarding)
    g : int, optional, default=None
       Statistical inefficiency to use if desired; if None, will be computed.

    TODO
    ----
    * Automatically determine 'ndiscard'.

    """
    # Get current dimensions.
    niterations = ncfile.variables['energies'].shape[0]
    nstates = ncfile.variables['energies'].shape[1]
    natoms = ncfile.variables['energies'].shape[2]

    # Extract energies.
    logger.info("Reading energies...")
    energies = ncfile.variables['energies']
    u_kln_replica = np.zeros([nstates, nstates, niterations], np.float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    logger.info("Done.")

    # Deconvolute replicas
    logger.info("Deconvoluting replicas...")
    u_kln = np.zeros([nstates, nstates, niterations], np.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    logger.info("Done.")

    # Compute total negative log probability over all iterations.
    u_n = np.zeros([niterations], np.float64)
    for iteration in range(niterations):
        u_n[iteration] = np.sum(np.diagonal(u_kln[:,:,iteration]))

    # Discard initial data to equilibration.
    u_kln_replica = u_kln_replica[:,:,ndiscard:]
    u_kln = u_kln[:,:,ndiscard:]
    u_n = u_n[ndiscard:]

    # Truncate to number of specified conforamtions to use
    if (nuse):
        u_kln_replica = u_kln_replica[:,:,0:nuse]
        u_kln = u_kln[:,:,0:nuse]
        u_n = u_n[0:nuse]

    # Subsample data to obtain uncorrelated samples
    N_k = np.zeros(nstates, np.int32)
    indices = timeseries.subsampleCorrelatedData(u_n, g=g) # indices of uncorrelated samples
    #print(u_n) # DEBUG
    #indices = range(0,u_n.size) # DEBUG - assume samples are uncorrelated
    N = len(indices) # number of uncorrelated samples
    N_k[:] = N
    u_kln = u_kln[:,:,indices]
    logger.info("number of uncorrelated samples:")
    logger.info(N_k)
    logger.info("")

    # Check for the expanded cutoff states, and subsamble as needed
    try:
        u_ln_full_raw = ncfile.variables['fully_interacting_expanded_cutoff_energies'][:].T #Its stored as nl, need in ln
        u_ln_non_raw = ncfile.variables['noninteracting_expanded_cutoff_energies'][:].T
        fully_interacting_u_ln = np.zeros(u_ln_full_raw.shape)
        noninteracting_u_ln = np.zeros(u_ln_non_raw.shape)
        # Deconvolute the fully interacting state
        for iteration in range(niterations):
            state_indices = ncfile.variables['states'][iteration,:]
            fully_interacting_u_ln[state_indices,iteration] = u_ln_full_raw[:,iteration]
            noninteracting_u_ln[state_indices,iteration] = u_ln_non_raw[:,iteration]
        # Discard non-equilibrated samples
        fully_interacting_u_ln = fully_interacting_u_ln[:,ndiscard:]
        fully_interacting_u_ln = fully_interacting_u_ln[:,indices]
        noninteracting_u_ln = noninteracting_u_ln[:,ndiscard:]
        noninteracting_u_ln = noninteracting_u_ln[:,indices]
        # Augment u_kln to accept the new state
        u_kln_new = np.zeros([nstates + 2, nstates + 2, N], np.float64)
        N_k_new = np.zeros(nstates + 2, np.int32)
        # Insert energies
        u_kln_new[1:-1,0,:] = fully_interacting_u_ln
        u_kln_new[1:-1,-1,:] = noninteracting_u_ln
        # Fill in other energies
        u_kln_new[1:-1,1:-1,:] = u_kln
        N_k_new[1:-1] = N_k
        # Notify users
        logger.info("Found expanded cutoff states in the energies!")
        logger.info("Free energies will be reported relative to them instead!")
        # Reset values, last step in case something went wrong so we dont overwrite u_kln on accident
        u_kln = u_kln_new
        N_k = N_k_new
    except:
        pass

    return u_kln, N_k, u_n


def initialize_MBAR(ncfile, u_kln=None, N_k=None):
    """
    Initialize MBAR for Free Energy and Enthalpy estimates, this may take a while.

    ncfile : NetCDF
       Input YANK netcdf file
    u_kln : array of numpy.float64, optional, default=None
       Reduced potential energies of the replicas; if None, will be extracted from the ncfile
    N_k : array of ints, optional, default=None
       Number of samples drawn from each kth replica; if None, will be extracted from the ncfile

    TODO
    ----
    * Ensure that the u_kln and N_k are decorrelated if not provided in this function

    """

    if u_kln is None or N_k is None:
        (u_kln, N_k, u_n) = extract_ncfile_energies(ncfile)

    # Initialize MBAR (computing free energy estimates, which may take a while)
    logger.info("Computing free energy differences...")
    mbar = MBAR(u_kln, N_k)

    return mbar


def estimate_free_energies(ncfile, mbar=None):
    """
    Estimate free energies of all alchemical states.

    Parameters
    ----------
    ncfile : NetCDF
       Input YANK netcdf file
    mbar : pymbar MBAR object, optional, default=None
       Initilized MBAR object from simulations; if None, it will be generated from the ncfile

    """

    # Create MBAR object if not provided
    if mbar is None:
        mbar = initialize_MBAR(ncfile)

    nstates = mbar.N_k.size

    # Get matrix of dimensionless free energy differences and uncertainty estimate.
    logger.info("Computing covariance matrix...")

    try:
        # pymbar 2
        (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences()
    except ValueError:
        # pymbar 3
        (Deltaf_ij, dDeltaf_ij, theta_ij) = mbar.getFreeEnergyDifferences()

    # Matrix of free energy differences
    logger.info("Deltaf_ij:")
    for i in range(nstates):
        str_row = ""
        for j in range(nstates):
            str_row += "%8.3f" % Deltaf_ij[i, j]
        logger.info(str_row)

    # Matrix of uncertainties in free energy difference (expectations standard
    # deviations of the estimator about the true free energy)
    logger.info("dDeltaf_ij:")
    for i in range(nstates):
        str_row = ""
        for j in range(nstates):
            str_row += "%8.3f" % dDeltaf_ij[i, j]
        logger.info(str_row)

    # Return free energy differences and an estimate of the covariance.
    return Deltaf_ij, dDeltaf_ij


def estimate_enthalpies(ncfile, mbar=None):
    """
    Estimate enthalpies of all alchemical states.

    Parameters
    ----------
    ncfile : NetCDF
       Input YANK netcdf file
    mbar : pymbar MBAR object, optional, default=None
       Initilized MBAR object from simulations; if None, it will be generated from the ncfile

    TODO
    ----
    * Check if there is an output/function name difference between pymbar 2 and 3
    """

    # Create MBAR object if not provided
    if mbar is None:
        mbar = initialize_MBAR(ncfile)

    nstates = mbar.N_k.size

    # Compute average enthalpies
    (f_k, df_k, H_k, dH_k, S_k, dS_k) = mbar.computeEntropyAndEnthalpy()

    return H_k, dH_k


def extract_u_n(ncfile):
    """
    Extract timeseries of u_n = - log q(X_n) from store file

    where q(X_n) = \pi_{k=1}^K u_{s_{nk}}(x_{nk})

    with X_n = [x_{n1}, ..., x_{nK}] is the current collection of replica configurations
    s_{nk} is the current state of replica k at iteration n
    u_k(x) is the kth reduced potential

    Parameters
    ----------
    ncfile : str
       The filename of the repex NetCDF file.

    Returns
    -------
    u_n : numpy array of numpy.float64
       u_n[n] is -log q(X_n)

    TODO
    ----
    Move this to repex.

    """

    # Get current dimensions.
    niterations = ncfile.variables['energies'].shape[0]
    nstates = ncfile.variables['energies'].shape[1]
    natoms = ncfile.variables['energies'].shape[2]

    # Extract energies.
    logger.info("Reading energies...")
    energies = ncfile.variables['energies']
    u_kln_replica = np.zeros([nstates, nstates, niterations], np.float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    logger.info("Done.")

    # Deconvolute replicas
    logger.info("Deconvoluting replicas...")
    u_kln = np.zeros([nstates, nstates, niterations], np.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    logger.info("Done.")

    # Compute total negative log probability over all iterations.
    u_n = np.zeros([niterations], np.float64)
    for iteration in range(niterations):
        u_n[iteration] = np.sum(np.diagonal(u_kln[:,:,iteration]))

    return u_n

# =============================================================================================
# SHOW STATUS OF STORE FILES
# =============================================================================================


def print_status(store_directory):
    """
    Print a quick summary of simulation progress.

    Parameters
    ----------
    store_directory : string
       The location of the NetCDF simulation output files.

    Returns
    -------
    success : bool
       True is returned on success; False if some files could not be read.

    """
    # Get NetCDF files
    phases = utils.find_phases_in_store_directory(store_directory)

    # Process each netcdf file.
    for phase, fullpath in phases.items():

        # Check that the file exists.
        if not os.path.exists(fullpath):
            # Report failure.
            logger.info("File %s not found." % fullpath)
            logger.info("Check to make sure the right directory was specified, and 'yank setup' has been run.")
            return False

        # Open NetCDF file for reading.
        logger.debug("Opening NetCDF trajectory file '%(fullpath)s' for reading..." % vars())
        ncfile = netcdf.Dataset(fullpath, 'r')

        # Read dimensions.
        niterations = ncfile.variables['positions'].shape[0]
        nstates = ncfile.variables['positions'].shape[1]
        natoms = ncfile.variables['positions'].shape[2]

        # Print summary.
        logger.info("%s" % phase)
        logger.info("  %8d iterations completed" % niterations)
        logger.info("  %8d alchemical states" % nstates)
        logger.info("  %8d atoms" % natoms)

        # TODO: Print average ns/day and estimated completion time.

        # Close file.
        ncfile.close()

    return True

# =============================================================================================
# ANALYZE STORE FILES
# =============================================================================================

def analyze(source_directory):
    """
    Analyze contents of store files to compute free energy differences.

    Parameters
    ----------
    source_directory : string
       The location of the NetCDF simulation storage files.

    Returns
    -------
    DeltaF, dDeltaF : float
       Dimensionless overall process
        free energy and statistical error

    """
    analysis_script_path = os.path.join(source_directory, 'analysis.yaml')
    if not os.path.isfile(analysis_script_path):
        err_msg = 'Cannot find analysis.yaml script in {}'.format(source_directory)
        logger.error(err_msg)
        raise RuntimeError(err_msg)
    with open(analysis_script_path, 'r') as f:
        analysis = yaml.load(f)
    phases = [phase_name for phase_name, sign in analysis]

    # Storage for different phases.
    data = dict()

    # Process each netcdf file.
    for phase in phases:
        ncfile_path = os.path.join(source_directory, phase + '.nc')

        # Open NetCDF file for reading.
        logger.info("Opening NetCDF trajectory file %(ncfile_path)s for reading..." % vars())
        try:
            ncfile = netcdf.Dataset(ncfile_path, 'r')

            logger.debug("dimensions:")
            for dimension_name in ncfile.dimensions.keys():
                logger.debug("%16s %8d" % (dimension_name, len(ncfile.dimensions[dimension_name])))

            # Read dimensions.
            niterations = ncfile.variables['positions'].shape[0]
            nstates = ncfile.variables['positions'].shape[1]
            logger.info("Read %(niterations)d iterations, %(nstates)d states" % vars())

            DeltaF_restraints = 0.0
            if 'metadata' in ncfile.groups:
                # Read phase direction and standard state correction free energy.
                # Yank sets correction to 0 if there are no restraints
                DeltaF_restraints = ncfile.groups['metadata'].variables['standard_state_correction'][0]

            # Choose number of samples to discard to equilibration
            MIN_ITERATIONS = 10 # minimum number of iterations to use automatic detection
            if niterations > MIN_ITERATIONS:
                from pymbar import timeseries
                u_n = extract_u_n(ncfile)
                u_n = u_n[1:] # discard initial frame of zero energies TODO: Get rid of initial frame of zero energies
                [nequil, g_t, Neff_max] = timeseries.detectEquilibration(u_n)
                nequil += 1 # account for initial frame of zero energies
                logger.info([nequil, Neff_max])
            else:
                nequil = 1  # discard first frame
                g_t = 1
                Neff_max = niterations

            # Examine acceptance probabilities.
            show_mixing_statistics(ncfile, cutoff=0.05, nequil=nequil)

            # Extract equilibrated, decorrelated energies, check for fully interacting state
            (u_kln, N_k, u_n) = extract_ncfile_energies(ncfile, ndiscard=nequil, g=g_t)

            # Create MBAR object to use for free energy and entropy states
            mbar = initialize_MBAR(ncfile, u_kln=u_kln, N_k=N_k)

            # Estimate free energies, use fully interacting state if present
            (Deltaf_ij, dDeltaf_ij) = estimate_free_energies(ncfile, mbar=mbar)

            # Estimate average enthalpies
            (DeltaH_i, dDeltaH_i) = estimate_enthalpies(ncfile, mbar=mbar)

            # Accumulate free energy differences
            entry = dict()
            entry['DeltaF'] = Deltaf_ij[0, -1]
            entry['dDeltaF'] = dDeltaf_ij[0, -1]
            entry['DeltaH'] = DeltaH_i[0, -1]
            entry['dDeltaH'] = dDeltaH_i[0, -1]
            entry['DeltaF_restraints'] = DeltaF_restraints
            data[phase] = entry

            # Get temperatures.
            ncvar = ncfile.groups['thermodynamic_states'].variables['temperatures']
            temperature = ncvar[0] * units.kelvin
            kT = kB * temperature

        finally:
            ncfile.close()

    # Compute free energy and enthalpy
    DeltaF = 0.0
    dDeltaF = 0.0
    DeltaH = 0.0
    dDeltaH = 0.0
    for phase, sign in analysis:
        DeltaF -= sign * (data[phase]['DeltaF'] + data[phase]['DeltaF_restraints'])
        dDeltaF += data[phase]['dDeltaF']**2
        DeltaH -= sign * (data[phase]['DeltaH'] + data[phase]['DeltaF_restraints'])
        dDeltaH += data[phase]['dDeltaH']**2
    dDeltaF = np.sqrt(dDeltaF)
    dDeltaH = np.sqrt(dDeltaH)

    # Attempt to guess type of calculation
    calculation_type = ''
    for phase in phases:
        if 'complex' in phase:
            calculation_type = ' of binding'
        elif 'solvent1' in phase:
            calculation_type = ' of solvation'

    # Print energies
    logger.info("")
    logger.info("Free energy{}: {:16.3f} +- {:.3f} kT ({:16.3f} +- {:.3f} kcal/mol)".format(
        calculation_type, DeltaF, dDeltaF, DeltaF * kT / units.kilocalories_per_mole,
        dDeltaF * kT / units.kilocalories_per_mole))
    logger.info("")

    for phase in phases:
        logger.info("DeltaG {:<25} : {:16.3f} +- {:.3f} kT".format(phase, data[phase]['DeltaF'],
                                                                   data[phase]['dDeltaF']))
        if data[phase]['DeltaF_restraints'] != 0.0:
            logger.info("DeltaG {:<25} : {:25.3f} kT".format('restraint',
                                                             data[phase]['DeltaF_restraints']))
    logger.info("")
    logger.info("Enthalpy{}: {:16.3f} +- {:.3f} kT ({:16.3f} +- {:.3f} kcal/mol)".format(
        calculation_type, DeltaH, dDeltaH, DeltaH * kT / units.kilocalories_per_mole,
        dDeltaH * kT / units.kilocalories_per_mole))

    return [DeltaF, dDeltaF]

# ==============================================================================
# Extract trajectory from NetCDF4 file
# ==============================================================================

def extract_trajectory(output_path, nc_path, state_index=None, replica_index=None,
                       start_frame=0, end_frame=-1, skip_frame=1, keep_solvent=True,
                       discard_equilibration=False, image_molecules=False):
    """Extract phase trajectory from the NetCDF4 file.

    Parameters
    ----------
    output_path : str
        Path to the trajectory file to be created. The extension of the file
        determines the format.
    nc_path : str
        Path to the NetCDF4 file containing the trajectory.
    state_index : int, optional
        The index of the alchemical state for which to extract the trajectory.
        One and only one between state_index and replica_index must be not None
        (default is None).
    replica_index : int, optional
        The index of the replica for which to extract the trajectory. One and
        only one between state_index and replica_index must be not None (default
        is None).
    start_frame : int, optional
        Index of the first frame to include in the trajectory (default is 0).
    end_frame : int, optional
        Index of the last frame to include in the trajectory. If negative, will
        count from the end (default is -1).
    skip_frame : int, optional
        Extract one frame every skip_frame (default is 1).
    keep_solvent : bool, optional
        If False, solvent molecules are ignored (default is True).
    discard_equilibration : bool, optional
        If True, initial equilibration frames are discarded (see the method
        pymbar.timeseries.detectEquilibration() for details, default is False).

    """
    # Check correct input
    if (state_index is None) == (replica_index is None):
        raise ValueError('One and only one between "state_index" and '
                         '"replica_index" must be specified.')
    if not os.path.isfile(nc_path):
        raise ValueError('Cannot find file {}'.format(nc_path))

    # Import simulation data
    try:
        nc_file = netcdf.Dataset(nc_path, 'r')

        # Extract topology and system serialization
        serialized_system = nc_file.groups['metadata'].variables['reference_system'][0]
        serialized_topology = nc_file.groups['metadata'].variables['topology'][0]

        # Determine if system is periodic
        from simtk import openmm
        reference_system = openmm.XmlSerializer.deserialize(str(serialized_system))
        is_periodic = reference_system.usesPeriodicBoundaryConditions()
        logger.info('Detected periodic boundary conditions: {}'.format(is_periodic))

        # Get dimensions
        n_iterations = nc_file.variables['positions'].shape[0]
        n_atoms = nc_file.variables['positions'].shape[2]
        logger.info('Number of iterations: {}, atoms: {}'.format(n_iterations, n_atoms))

        # Determine frames to extract
        if start_frame <= 0:
            # TODO yank saves first frame with 0 energy!
            start_frame = 1
        if end_frame < 0:
            end_frame = n_iterations + end_frame + 1
        frame_indices = range(start_frame, end_frame, skip_frame)
        if len(frame_indices) == 0:
            raise ValueError('No frames selected')
        logger.info('Extracting frames from {} to {} every {}'.format(
            start_frame, end_frame, skip_frame))

        # Discard equilibration samples
        if discard_equilibration:
            u_n = extract_u_n(nc_file)[frame_indices]
            n_equil, g, n_eff = timeseries.detectEquilibration(u_n)
            logger.info(("Discarding initial {} equilibration samples (leaving {} "
                         "effectively uncorrelated samples)...").format(n_equil, n_eff))
            frame_indices = frame_indices[n_equil:-1]

        # Extract state positions and box vectors
        positions = np.zeros((len(frame_indices), n_atoms, 3))
        if is_periodic:
            box_vectors = np.zeros((len(frame_indices), 3, 3))
        if state_index is not None:
            logger.info('Extracting positions of state {}...'.format(state_index))

            # Deconvolute state indices
            state_indices = np.zeros(len(frame_indices))
            for i, iteration in enumerate(frame_indices):
                replica_indices = nc_file.variables['states'][iteration, :]
                state_indices[i] = np.where(replica_indices == state_index)[0][0]

            # Extract state positions and box vectors
            for i, iteration in enumerate(frame_indices):
                replica_index = state_indices[i]
                positions[i, :, :] = nc_file.variables['positions'][iteration, replica_index, :, :].astype(np.float32)
                if is_periodic:
                    box_vectors[i, :, :] = nc_file.variables['box_vectors'][iteration, replica_index, :, :].astype(np.float32)

        else:  # Extract replica positions and box vectors
            logger.info('Extracting positions of replica {}...'.format(replica_index))

            for i, iteration in enumerate(frame_indices):
                positions[i, :, :] = nc_file.variables['positions'][iteration, replica_index, :, :].astype(np.float32)
                if is_periodic:
                    box_vectors[i, :, :] = nc_file.variables['box_vectors'][iteration, replica_index, :, :].astype(np.float32)
    finally:
        nc_file.close()

    # Create trajectory object
    logger.info('Creating trajectory object...')
    topology = utils.deserialize_topology(serialized_topology)
    trajectory = mdtraj.Trajectory(positions, topology)
    if is_periodic:
        trajectory.unitcell_vectors = box_vectors

    # Force periodic boundary conditions to molecules positions
    if image_molecules:
        logger.info('Applying periodic boundary conditions to molecules positions...')
        trajectory.image_molecules(inplace=True)

    # Remove solvent
    if not keep_solvent:
        logger.info('Removing solvent molecules...')
        trajectory = trajectory.remove_solvent()

    # Detect format
    extension = os.path.splitext(output_path)[1][1:]  # remove dot
    try:
        save_function = getattr(trajectory, 'save_' + extension)
    except AttributeError:
        raise ValueError('Cannot detect format from extension of file {}'.format(output_path))

    # Create output directory and save trajectory
    logger.info('Creating trajectory file: {}'.format(output_path))
    output_dir = os.path.dirname(output_path)
    if output_dir != '' and not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    save_function(output_path)
