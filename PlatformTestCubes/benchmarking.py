from __future__ import print_function
"""
Modified variant of benchmark.py from https://github.com/pandegroup/openmm/blob/master/examples/benchmark.py
"""
import os
import sys
from datetime import datetime

import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as unit


def timeIntegration(context, steps, initialSteps):
    """Integrate a Context for a specified number of steps, then return how many seconds it took."""
    context.getIntegrator().step(initialSteps)  # Make sure everything is fully initialized
    context.getState(getEnergy=True)
    start = datetime.now()
    context.getIntegrator().step(steps)
    context.getState(getEnergy=True)
    end = datetime.now()
    elapsed = end - start
    return elapsed.seconds + elapsed.microseconds * 1e-6


def run_platform_benchmarks(options, stream=None):
    if stream is None:
        stream = sys.stdout
    platforms = [mm.Platform.getPlatform(i) for i in range(mm.Platform.getNumPlatforms())]
    for platform in platforms:
        stream.write("Performing Benchmarks on {}\n".format(platform.getName()))
        options.__dict__.update({"platform": platform})
        for test_type in ["rf", "pme", "amoebagk", "amoebapme"]:
            stream.write("Running Benchmarks for {} using {}\n".format(test_type, platform.getName()))
            run_benchmark("{}_{}".format(platform.getName(), test_type), options, stream=stream)
            stream.write("\n\n")  # add line breaks to seperate runs
            stream.flush()


def run_benchmark(test_name, options, stream=None):
    """Perform a single benchmarking simulation."""
    if stream is None:
        stream = sys.stdout
    benchmark_data_dir = os.path.join(os.path.dirname(__file__), "data")

    explicit = any([x in test_name for x in ["rf", "pme", "amoebapme"]])
    amoeba = any([x in test_name for x in ["amoebagk", "amoebapme"]])
    hydrogenMass = None
    if amoeba:
        stream.write('Test: %s (epsilon=%g)\n' % (test_name, options.amoeba_target_epsilon))
    elif test_name == 'pme':
        stream.write('Test: pme (cutoff=%g)\n' % options.cutoff)
    else:
        stream.write('Test: %s\n' % test_name)

    # Create the System.

    if amoeba:
        constraints = None
        epsilon = float(options.amoeba_target_epsilon)
        if explicit:
            ff = app.ForceField('amoeba2009.xml')
            pdb = app.PDBFile(os.path.join(benchmark_data_dir, '5dfr_solv-cube_equil.pdb'))
            cutoff = 0.7 * unit.nanometers
            vdwCutoff = 0.9 * unit.nanometers
            system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, vdwCutoff=vdwCutoff, constraints=constraints, ewaldErrorTolerance=0.00075, mutualInducedTargetEpsilon=epsilon, polarization=options.polarization)
        else:
            ff = app.ForceField(
                'amoeba2009.xml',
                'amoeba2009_gk.xml'
            )
            pdb = app.PDBFile(os.path.join(benchmark_data_dir, '5dfr_minimized.pdb'))
            cutoff = 2.0 * unit.nanometers
            vdwCutoff = 1.2 * unit.nanometers
            system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=constraints, mutualInducedTargetEpsilon=epsilon, polarization=options.polarization)
        for f in system.getForces():
            if isinstance(f, mm.AmoebaMultipoleForce) or isinstance(f, mm.AmoebaVdwForce) or isinstance(f, mm.AmoebaGeneralizedKirkwoodForce) or isinstance(f, mm.AmoebaWcaDispersionForce):
                f.setForceGroup(1)
        dt = 0.002 * unit.picoseconds
        integ = mm.MTSIntegrator(dt, [(0, 2), (1, 1)])
    else:
        if explicit:
            ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
            pdb = app.PDBFile(os.path.join(benchmark_data_dir, '5dfr_solv-cube_equil.pdb'))
            if 'pme' in test_name:
                method = app.PME
                cutoff = options.cutoff
            else:
                method = app.CutoffPeriodic
                cutoff = 1 * unit.nanometers
            friction = 1 * (1 / unit.picoseconds)
        else:
            ff = app.ForceField(
                'amber99sb.xml',
                'amber99_obc.xml'
            )
            pdb = app.PDBFile(os.path.join(benchmark_data_dir, '5dfr_minimized.pdb'))
            method = app.CutoffNonPeriodic
            cutoff = 2 * unit.nanometers
            friction = 91 * (1 / unit.picoseconds)
        if options.use_heavy_hydrogens:
            dt = 0.005 * unit.picoseconds
            constraints = app.AllBonds
            hydrogenMass = 4 * unit.amu
        else:
            dt = 0.002 * unit.picoseconds
            constraints = app.HBonds
            hydrogenMass = None
        system = ff.createSystem(pdb.topology, nonbondedMethod=method, nonbondedCutoff=cutoff, constraints=constraints, hydrogenMass=hydrogenMass)
        integ = mm.LangevinIntegrator(300 * unit.kelvin, friction, dt)
    stream.write('Step Size: %g fs\n' % dt.value_in_unit(unit.femtoseconds))
    properties = {}
    initialSteps = 5
    if options.precision is not None and options.platform.getName() in ('CUDA', 'OpenCL'):
        properties['Precision'] = options.precision

    # Run the simulation.

    integ.setConstraintTolerance(1e-5)
    try:
        if len(properties) > 0:
            context = mm.Context(system, integ, options.platform, properties)
        else:
            context = mm.Context(system, integ, options.platform)
    except Exception as e:
        stream.write("Unable to Benchmark '{}'\n".format(test_name))
        stream.write("ERROR: {}\n".format(str(e)))
        return
    context.setPositions(pdb.positions)
    context.setVelocitiesToTemperature(300 * unit.kelvin)
    steps = 20
    while True:
        time = timeIntegration(context, steps, initialSteps)
        if time >= 0.5 * options.seconds:
            break
        if time < 0.5:
            steps = int(steps * 1.0 / time)  # Integrate enough steps to get a reasonable estimate for how many we'll need.
        else:
            steps = int(steps * options.seconds / time)
    stream.write('Integrated %d steps in %g seconds\n' % (steps, time))
    stream.write('%g ns/day\n' % (dt * steps * 86400 / time).value_in_unit(unit.nanoseconds))
