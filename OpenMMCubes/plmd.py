#############################################################################
# Copyright (C) 2017 OpenEye Scientific Software, Inc.
#############################################################################

import io, os, time, logging, traceback
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


def MakeAtomMaskLigNonH( plMask ):
    """return an integer list the size of plMask, making >0 only the ligand nonH.
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
        # otherwise, no restraint
        else:
            atomMask.append(0)
    return atomMask


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


def MakeAtomMaskLigKeyCA( plMask ):
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
def RestrMin(topology, system, positions, PLMask, restrwt=2.0, steps=100):
    """" Warms up the systems to the secified temp, restraining nonWat nonH.
           topology: the openmm Topology object for the whole system in the periodic box
           system: the openmm System object for the whole system in the periodic box
           PLmask: the protein-ligand mask used to assign restraints
           restrwt: the weight in kcal/mol/ang^2 to be applied to restrained atoms
           runTime: the length in picoseconds of the warmup run
           targetTemp: the target temperature in Kelvin
           """
    print('RestrMin: Minimization for %d steps with %3.1f kcal/mol/ang^2 restraints on all non-water non-Hydrogens'
          % (steps, restrwt))
    overall_timer = LoggingStopwatch()
    stage_timer = LoggingStopwatch()

    stage_timer.TimeCheck('RestrMin: building system and simulation objects')
    restrMask = MakeAtomMaskNonSolventNonH( PLMask)
    restrStrongNonSolvNonH = MakeOpenMMRestraintForceObj( positions, restrMask, restrwt)
    system.addForce( restrStrongNonSolvNonH)

    integrator = openmm.LangevinIntegrator(296 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
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


if __name__ == "__main__":
    sys.exit(main(sys.argv))
