#############################################################################
# Copyright (C) 2017 OpenEye Scientific Software, Inc.
#############################################################################

import io, os, time, traceback
from openeye import oechem
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app


#############################################################################
# make a atom mask file for protein-ligand MD based on a molecule for the solvated complex
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


if __name__ == "__main__":
    sys.exit(main(sys.argv))
