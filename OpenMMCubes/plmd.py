#############################################################################
# Copyright (C) 2017 OpenEye Scientific Software, Inc.
#############################################################################

import io, os, time, logging, traceback
from sys import stdout
from openeye import oechem
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app

FORMAT = "%(asctime)-15s %(epoch).2f | %(message)s"
TIMEFORMAT = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(datefmt=TIMEFORMAT, format=FORMAT, level=logging.INFO)

class LoggingStopwatch:
    def __init__(self):
        self.start_time = time.time()

    def Start(self):
        self.start_time = time.time()

    def Stop(self, stagename):
        end_time = time.time()
        duration = end_time - self.start_time
        message  = "%.2f CPU sec: " % duration
        logging.info(message + stagename, extra={"epoch":end_time})

    def TimeCheck(self, stagename):
        end_time = time.time()
        duration = end_time - self.start_time
        message  = "%.2f CPU sec: " % duration
        logging.info(message + stagename, extra={"epoch":end_time})
        self.start_time = time.time()



def MakeOpenMMRestraintForceObj( particlePositions, AtomMask, restWt=0.5 ):
    """return an openmm CustomExternalForce object based on passed xyz coords, restraining
       the atoms specified >0 in the AtomMask.
       particlePositions: the openmm xyz coordinate object for the whole system
       AtomMask: a list of integers, 0 for unrestrained atoms and positive for restrained atoms.
       restWt: restraint weight in kcal/mol/ang^2 (default 0.05)"""
    # define the custom force to restrain atoms to their starting positions
    force_restr = openmm.CustomExternalForce('k_restr*( (x-x0)^2 + (y-y0)^2 + (z-z0)^2 )')
    # Add the restraint weight as a global parameter in kcal/mol/A^2
    force_restr.addGlobalParameter("k_restr", restWt*unit.kilocalories_per_mole/unit.angstroms**2)
    # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
    force_restr.addPerParticleParameter("x0")
    force_restr.addPerParticleParameter("y0")
    force_restr.addPerParticleParameter("z0")
    # loop over all atoms (particles), zipping the mask list with the atoms
    for i, (xyz_tuple, irestr) in enumerate(zip(particlePositions, AtomMask)):
    # add the restraint if the AtomMask variable irestr is >0
        if irestr > 0:
            force_restr.addParticle(i, xyz_tuple)
    return force_restr


def MakeAtomMaskNonSolventNonH( plMask ):
    """return an integer list the size of plMask, making >0 only non-Solvent nonH.
       plMask: a dictionary of three one-per-atom lists: ['PLMask'], where each atom is
          designated "Protein", "Water", "Ion", or "Other; ['AtNum'] and ['AtName'],
          which are the atom's atomic number and pdb atom name, respectively."""
    # loop over all atoms (particles), zipping the mask list with the atoms
    atomMask = []
    for i, (atype, atnum, atname, asite) in enumerate(
            zip(plMask['PLMask'], plMask['AtNum'], plMask['AtName'], plMask['ActSite'])):
    # add the restraint if the atom is non-Water and nonH
        if atype != "Water" and atype != "Ion" and atnum > 1:
            atomMask.append(1)
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask


def MakeAtomMaskProteinNonH( plMask ):
    """return an integer list the size of plMask, making >0 only the protein nonH.
       plMask: a dictionary of three one-per-atom lists: ['PLMask'], where each atom is
          designated "Protein", "Water", or "Other; ['AtNum'] and ['AtName'], which are
          the atom's atomic number and pdb atom name, respectively."""
    # loop over all atoms (particles), zipping the mask list with the atoms
    atomMask = []
    for i, (atype, atnum, atname, asite) in enumerate(
            zip(plMask['PLMask'], plMask['AtNum'], plMask['AtName'], plMask['ActSite'])):
    # add the restraint if the atom is Protein and nonH
        if atype == "Protein" and atnum > 1:
            atomMask.append(1)
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask


def MakeAtomMaskProteinCAlpha( plMask ):
    """return an integer list the size of plMask, making >0 only the protein alpha carbons.
       plMask: a dictionary of three one-per-atom lists: ['PLMask'], where each atom is
          designated "Protein", "Water", or "Other; ['AtNum'] and ['AtName'], which are
          the atom's atomic number and pdb atom name, respectively."""
    # loop over all atoms (particles), zipping the mask list with the atoms
    atomMask = []
    for i, (atype, atnum, atname, asite) in enumerate(
            zip(plMask['PLMask'], plMask['AtNum'], plMask['AtName'], plMask['ActSite'])):
    # add the restraint if the atom is protein and atom name is CA
        if atype == "Protein" and atname == "CA":
            atomMask.append(1)
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask


def MakeAtomMaskLigandNonH( plMask ):
    """return an integer list the size of plMask, making >0 only the ligand nonH.
       plMask: a dictionary of three one-per-atom lists: ['PLMask'], where each atom is
          designated "Protein", "Water", or "Other; ['AtNum'] and ['AtName'], which are
          the atom's atomic number and pdb atom name, respectively."""
    # loop over all atoms (particles), zipping the mask list with the atoms
    atomMask = []
    for i, (atype, atnum, atname, asite) in enumerate(
            zip(plMask['PLMask'], plMask['AtNum'], plMask['AtName'], plMask['ActSite'])):
    # add the restraint if the atom is Other and nonH
        if atype == "Other" and atnum > 1:
            atomMask.append(1)
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask


def MakeAtomMaskLigandKeyCA( plMask ):
    """return an integer list the size of plMask, making >0
       only the ligand nonH and specified (in the plMask) alpha-carbons in the protein.
       plMask: a dictionary of three one-per-atom lists: ['PLMask'], where each atom is
          designated "Protein", "Water", or "Other; ['AtNum'] and ['AtName'], which are
          the atom's atomic number and pdb atom name, respectively."""
    # loop over all atoms (particles), zipping the mask list with the atoms
    atomMask = []
    for i, (atype, atnum, atname, asite) in enumerate(
            zip(plMask['PLMask'], plMask['AtNum'], plMask['AtName'], plMask['ActSite'])):
    # add the restraint if the atom is non-Water and nonH
        if atype == "Other" and atnum > 1:
            atomMask.append(1)
        # add the restraint if the atom is an active-site CA
        elif atname == "CA" and asite > 0:
            atomMask.append(1)
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask


def MakeAtomMaskKeyCA( plMask ):
    """return an integer list the size of plMask, making >0
       only the specified (in the plMask) alpha-carbons in the protein.
       plMask: a dictionary of three one-per-atom lists: ['PLMask'], where each atom is
          designated "Protein", "Water", or "Other; ['AtNum'] and ['AtName'], which are
          the atom's atomic number and pdb atom name, respectively."""
    # loop over all atoms (particles), zipping the mask list with the atoms
    atomMask = []
    for i, (atype, atnum, atname, asite) in enumerate(
            zip(plMask['PLMask'], plMask['AtNum'], plMask['AtName'], plMask['ActSite'])):
    # add the restraint if the atom is an active-site CA
        if atname == "CA" and asite > 0:
            atomMask.append(1)
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask

#############################################################################
# make an atom mask file for protein-ligand MD based on a molecule for the solvated complex
#############################################################################
def protLigMask(mol, actSiteResNumTag):
    if oechem.OEHasSDData( mol, actSiteResNumTag):
        actSiteResNumsStr = oechem.OEGetSDData(mol, actSiteResNumTag)
        actSiteResNums = [ int(resNum) for resNum in actSiteResNumsStr.split()]
    else:
        actSiteResNums = []
    atomProps = { 'PLMask':[], 'AtNum':[], 'ResNum':[], 'AtName':[], 'ActSite':[] }
    for atom in mol.GetAtoms():
        atomProps['AtNum'].append( atom.GetAtomicNum() )
        atomProps['AtName'].append( atom.GetName().strip() )
        res= oechem.OEAtomGetResidue(atom)
        resname= res.GetName()
        resindx= oechem.OEGetResidueIndex(res)
        atype= "Other"
        if oechem.OEIsStandardProteinResidue( resindx):
            atype= "Protein"
        elif resname=="ACE" or resname=="NME":
            atype= "Protein"
        elif resindx== oechem.OEResidueIndex_HOH:
            atype= "Water"
        elif atom.GetDegree()<1:
            atype= "Ion"
        atomProps['PLMask'].append( atype )
        resnum= res.GetResidueNumber()
        atomProps['ResNum'].append( resnum )
        if resnum in actSiteResNums:
            isActSiteRes= 1
        else:
            isActSiteRes= 0
        atomProps['ActSite'].append( isActSiteRes)
    return atomProps

#############################################################################
# minimize a solvated protein-ligand complex with restraints on non-solvent non-hydrogens
#############################################################################
def RestrMin( openmmStuff, options):
    """" Warms up the systems to the secified temp, restraining nonWat nonH.

         input parameters:
           openmmStuff: dict containing the OpenMM topology, system, state,
           initial positions, protein-ligand mask for restraints, and an idtag string.

           options: dict containing minimization parameters such as restraint weight
           and number of minimization steps.
           """
    overall_timer = LoggingStopwatch()
    stage_timer = LoggingStopwatch()

    topology = openmmStuff['topology']
    system = openmmStuff['system']
    positions = openmmStuff['positions']
    PLMask = openmmStuff['PLmask']
    restrwt = options['restraintWt']
    steps = options['steps']

    print('RestrMin: Minimization for %d steps with %3.1f kcal/mol/ang^2 restraints on all non-water non-Hydrogens'
          % (steps, restrwt))

    restrMask = MakeAtomMaskNonSolventNonH( PLMask)
    restrStrongNonSolvNonH = MakeOpenMMRestraintForceObj( positions, restrMask, restrwt)
    system.addForce( restrStrongNonSolvNonH)

    integrator = openmm.LangevinIntegrator( 300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    # minimize for minSteps, getting before-and-after energies
    state = simulation.context.getState(getEnergy=True)
    stage_timer.TimeCheck('RestrMin: Initial energy %.4f kcal/mol'
          % state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
    stage_timer.TimeCheck('RestrMin: Minimizing for %d steps' % steps)
    simulation.minimizeEnergy( maxIterations=steps)
    state = simulation.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
    stage_timer.TimeCheck('RestrMin: Final energy %.4f kcal/mol'
          % state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))

    return state

#############################################################################
# warm up a solvated protein-ligand complex with restraints on non-solvent non-hydrogens
#############################################################################
def RestrWarmupNVT( openmmStuff, options):
    """" Warms up the systems to the secified temp, restraining nonWat nonH.

         input parameters:
           openmmStuff: dict containing the OpenMM topology, system, state,
           initial positions, protein-ligand mask for restraints, and an idtag string.

           options: dict containing minimization parameters such as restraint weight,
           length of MD run in picoseconds, and target temperature.
           """
    topology = openmmStuff['topology']
    system = openmmStuff['system']
    PLMask = openmmStuff['PLmask']
    if openmmStuff['state']:
        state = openmmStuff['state']
        positions = state.getPositions()
    else:
        positions = openmmStuff['positions']

    picosec = options['picosec']
    temperature = options['temperature']
    restraintWt = options['restraintWt']

    print('RestrWarmupNVT: Warm up for %d picoseconds with %3.1f kcal/mol/ang^2 restraints on all non-water non-Hydrogens'
          % (picosec, restraintWt))
    overall_timer = LoggingStopwatch()
    stage_timer = LoggingStopwatch()

    print('RestrWarmupNVT: building system and simulation objects')
    restrMask = MakeAtomMaskNonSolventNonH( PLMask)
    restrStrongNonSolvNonH = MakeOpenMMRestraintForceObj( positions, restrMask, restraintWt)
    system.addForce( restrStrongNonSolvNonH)

    stepLen = 0.002
    tempFricCoeff= 0.5
    integrator = openmm.LangevinIntegrator(temperature*unit.kelvin,
                                           tempFricCoeff/unit.picosecond,
                                           stepLen*unit.picoseconds)
    simulation = app.Simulation(topology, system, integrator)

    # set positions and update the periodic box vectors
    simulation.context.setPositions( positions)
    if openmmStuff['state']:
        xvec, yvec, zvec = state.getPeriodicBoxVectors()
        simulation.context.setPeriodicBoxVectors( xvec, yvec, zvec)
    stage_timer.TimeCheck("RestrWarmupNVT: system and simulation objects created and initialized")

    # set up and run dynamics
    print('RestrWarmupNVT: Computations will be done on platform: ' + simulation.context.getPlatform().getName() )
    reportFreq = 100
    outfname = 'output/'+openmmStuff['idtag']+'-warmup.log'
    fileReporter = app.StateDataReporter( outfname, reportFreq, step=True,
            time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True)
    stdoutReporter = app.StateDataReporter( stdout, reportFreq, step=True,
            time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True)
    simulation.reporters.append( fileReporter)
    simulation.reporters.append( stdoutReporter)

    totalTimeSteps= picosec/stepLen
    print('RestrWarmupNVT: starting MD for %4.3f picoseconds (%d timesteps)' % (picosec, totalTimeSteps))
    simulation.step(totalTimeSteps)
    stage_timer.TimeCheck("RestrWarmupNVT: finishing MD run")

    state= simulation.context.getState( getPositions=True, getEnergy=True,
                                       getVelocities=True, enforcePeriodicBox=True )
    print('RestrWarmupNVT: system energy  after MD: %s'
          % state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)*unit.kilocalories_per_mole)

    return state

#############################################################################
# equilibrate a solvated protein-ligand complex with restraints
#############################################################################
def RestrEquil( openmmStuff, options):
    """" Equilibrates the NPT system at the specified temp, restraining nonWat nonH.

         input parameters:
           openmmStuff: dict containing the OpenMM topology, system, state,
           initial positions, protein-ligand mask for restraints, and an idtag string.

           options: dict containing minimization parameters such as restraint weight,
           length of MD run in picoseconds, and target temperature.
           """
    topology = openmmStuff['topology']
    system = openmmStuff['system']
    PLMask = openmmStuff['PLmask']
    state = openmmStuff['state']

    picosec = options['picosec']
    temperature = options['temperature']
    restraintWt = options['restraintWt']
    maskType = options['restraintType']
    label = options['label']

    print('RestrEquil: Equilibrate for %d picoseconds with %3.1f kcal/mol/ang^2 restraints on all %s'
          % (picosec, restraintWt, maskType))
    overall_timer = LoggingStopwatch()
    stage_timer = LoggingStopwatch()

    print('RestrEquil: building system and simulation objects')
    if maskType=='NonSolventNonH':
        restrMask = MakeAtomMaskNonSolventNonH( PLMask)
    if maskType=='ProteinNonH':
        restrMask = MakeAtomMaskProteinNonH( PLMask)
    if maskType=='ProteinCAlpha':
        restrMask = MakeAtomMaskProteinCAlpha( PLMask)
    if maskType=='LigandNonH':
        restrMask = MakeAtomMaskLigandNonH( PLMask)
    else:
        restrMask = None
    if restrMask and restraintWt>0.0:
        restrStrongNonSolvNonH = MakeOpenMMRestraintForceObj( state.getPositions(), restrMask, restraintWt)
        system.addForce( restrStrongNonSolvNonH)

    # make this an NPT (constant pressure) simulation by adding a barostat
    system.addForce( openmm.MonteCarloBarostat( 1*unit.bar, temperature*unit.kelvin))

    stepLen = 0.002
    tempFricCoeff= 1.0
    integrator = openmm.LangevinIntegrator(temperature*unit.kelvin,
                                           tempFricCoeff/unit.picosecond,
                                           stepLen*unit.picoseconds)
    simulation = app.Simulation(topology, system, integrator)

    # set the state to the saved state.
    simulation.context.setState( state)
    stage_timer.TimeCheck("RestrEquil: system and simulation objects created and initialized")

    # set up and run dynamics
    print('RestrEquil: Computations will be done on platform: ' + simulation.context.getPlatform().getName() )
    outfname = 'output/'+openmmStuff['idtag']+label
    # Determine reporting frequency in steps and optionally snapshot frequency and reporter
    snapFreq = options['snapFreq']
    if snapFreq > 0:
        print( 'requested snapshot frequency in ps is', snapFreq)
        reportFreq= int( snapFreq/stepLen)
        simulation.reporters.append(app.DCDReporter(outfname+'.dcd', reportFreq))
    else:
        reportFreq = int(1.0/stepLen)
    # set reporters
    print( 'reporting frequency is', reportFreq)
    simulation.reporters.append(app.StateDataReporter(outfname+'.log', reportFreq, step=True, time=True,
          potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
          temperature=True, volume=False, density=True ))
    simulation.reporters.append(app.StateDataReporter(stdout, reportFreq, step=True, time=True,
          potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
          temperature=True, volume=False, density=True ))

    totalTimeSteps= int(picosec/stepLen)
    print( 'totalTimeSteps is', totalTimeSteps)
    print('RestrEquil: starting MD for %4.3f picoseconds (%d timesteps)' % (picosec, totalTimeSteps))
    simulation.step(totalTimeSteps)
    stage_timer.TimeCheck("RestrEquil: finishing MD run")

    state= simulation.context.getState( getPositions=True, getEnergy=True,
                                       getVelocities= True, enforcePeriodicBox= True )
    print('RestrEquil: system energy  after MD: %s'
          % state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)*unit.kilocalories_per_mole)

    return state


if __name__ == "__main__":
    sys.exit(main(sys.argv))
