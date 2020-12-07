################################################################
#  Copyright (C) 2015 OpenEye Scientific Software, Inc.
################################################################
# import os, sys
# import numpy as np

from openeye import oechem

def GetInteractionString(inter):

    fragstrs = []
    for frag in [inter.GetBgnFragment(), inter.GetEndFragment()]:
        if frag is None:
            continue
        fragtype = frag.GetComponentType()
        if fragtype == oechem.OELigandInteractionHintComponent():
            fragstrs.append("ligand:" + " ".join(sorted(str(a) for a in frag.GetAtoms())))
        if fragtype == oechem.OEProteinInteractionHintComponent():
            fragstrs.append("protein:" + " ".join(sorted(str(a) for a in frag.GetAtoms())))

    return " ".join(sorted(f for f in fragstrs))


def PrintResidueInteractions(asite, residue):

    ligatomnames = set()
    for inter in asite.GetInteractions(oechem.OEHasResidueInteractionHint(residue)):
        ligfrag = inter.GetFragment(oechem.OELigandInteractionHintComponent())
        if ligfrag is None:
            continue
        for latom in ligfrag.GetAtoms():
            ligatomnames.add(str(latom))

    if len(ligatomnames) != 0:
        print(GetResidueName(residue), ": ", " ".join(sorted(a for a in ligatomnames)))


def PrintLigandAtomInteractions(asite, atom):

    resnames = set()
    for inter in asite.GetInteractions(oechem.OEHasInteractionHint(atom)):
        profrag = inter.GetFragment(oechem.OEProteinInteractionHintComponent())
        if profrag is None:
            continue
        for patom in profrag.GetAtoms():
            residue = oechem.OEAtomGetResidue(patom)
            resnames.add(GetResidueName(residue))

    if len(resnames) != 0:
        print(atom, ":", " ".join(sorted(r for r in resnames)))


def WriteInteractionHintsToString(hint_object,protein_name,ligand_name):
    hint_str = 'Protein-ligand interactions for {} {}\n'.format(protein.GetTitle(), ligand.GetTitle())
    hint_str += "Number of interactions: {}\n".format(asite.NumInteractions())
    hints = []
    for itype in oechem.OEGetActiveSiteInteractionHintTypes():
        numinters = asite.NumInteractions(oechem.OEHasInteractionHintType(itype))
        if numinters == 0:
            continue
        hint_str += "%d %s :\n" % (numinters, itype.GetName())

        inters = [s for s in asite.GetInteractions(oechem.OEHasInteractionHintType(itype))]
        hints.append(inters)
        hint_str += "\n".join(sorted(GetInteractionString(s) for s in inters))
        hint_str += "\n"

    return hint_str


def DictAllowedInteractionStrings(hint_object,allowedInteractionStrings):
    hint_intnstr = dict()
    for itype in oechem.OEGetActiveSiteInteractionHintTypes():
        numinters = hint_object.NumInteractions(oechem.OEHasInteractionHintType(itype))
        hname = itype.GetName()
        if numinters == 0 or hname not in allowedInteractionStrings:
            continue
        inters = [s for s in hint_object.GetInteractions(oechem.OEHasInteractionHintType(itype))]
        hint_intnstr[hname] = [GetInteractionString(s) for s in inters]
    return hint_intnstr


def AddFrameHintOverlapToRefHints(refHintDictDict,frameHintDict,frameId):
    frame_hintKeys = set(frameHintDict.keys())
    for key in refHintDictDict.keys():
        #print('refHintDictDict key:',key)
        # Pass over refHintDictDict keys not present in frameHintDict keys
        if key not in frame_hintKeys:
            continue
        # The key is in both refHintDictDict and frameHintDict so
        # append refHintDictDict vector for interactions common to both
        frame_inters = frameHintDict[key]
        #print(frame_inters)
        for inter in refHintDictDict[key].keys():
            if inter in frame_inters:
                #print('  interaction present in both:', inter)
                refHintDictDict[key][inter] += 1
                #print('    hint_traj',key,inter, hint_traj[key][inter])
    return True


allIntnTypes = [
'bio:active-site:contact',
'bio:active-site:clash',
'bio:active-site:covalent',
'bio:active-site:hbond:protein2ligand',
'bio:active-site:hbond:protein2ligand-charged',
'bio:active-site:hbond:ligand2protein',
'bio:active-site:hbond:ligand2protein-charged',
'bio:active-site:hbond:clash-acc-acc',
'bio:active-site:hbond:clash-acc-acc-charged',
'bio:active-site:hbond:clash-don-don',
'bio:active-site:hbond:clash-don-don-charged',
'bio:active-site:hbond:non-ideal-protein2ligand',
'bio:active-site:hbond:non-ideal-protein2ligand-charged',
'bio:active-site:hbond:non-ideal-ligand2protein',
'bio:active-site:hbond:non-ideal-ligand2protein-charged',
'bio:active-site:hbond:ligand-intra-molecular',
'bio:active-site:hbond:ligand-intra-molecular-charged',
'bio:active-site:hbond:protein-intra-molecular',
'bio:active-site:hbond:protein-intra-molecular-charged',
'bio:active-site:hbond:unpaired-ligand-donor',
'bio:active-site:hbond:unpaired-ligand-donor-charged',
'bio:active-site:hbond:unpaired-ligand-acceptor',
'bio:active-site:hbond:unpaired-ligand-acceptor-charged',
'bio:active-site:hbond:unpaired-protein-donor',
'bio:active-site:hbond:unpaired-protein-donor-charged',
'bio:active-site:hbond:unpaired-protein-acceptor',
'bio:active-site:hbond:unpaired-protein-acceptor-charged',
'bio:active-site:chelator:ligand-chelates',
'bio:active-site:chelator:protein-chelates',
'bio:active-site:chelator:ligand-intra-molecular',
'bio:active-site:chelator:protein-intra-molecular',
'bio:active-site:chelator:ligand-intra-molecular',
'bio:active-site:chelator:protein-intra-molecular',
'bio:active-site:salt-bridge:ligand-protein+',
'bio:active-site:salt-bridge:ligand+protein-',
'bio:active-site:salt-bridge:clash',
'bio:active-site:salt-bridge:unpaired-ligand+',
'bio:active-site:salt-bridge:unpaired-ligand-',
'bio:active-site:salt-bridge:unpaired-protein+',
'bio:active-site:salt-bridge:unpaired-protein-',
'bio:active-site:stacking:t',
'bio:active-site:stacking:pi',
'bio:active-site:cationpi:ligandpi',
'bio:active-site:cationpi:proteinpi',
'bio:active-site:halogen:ligand-electrophile',
'bio:active-site:halogen:ligand-nucleophile',
'bio:active-site:halogen:protein-electrophile',
'bio:active-site:halogen:protein-nucleophile',
]
goodIntnTypes = [
'bio:active-site:contact',
'bio:active-site:covalent',
'bio:active-site:hbond:protein2ligand',
'bio:active-site:hbond:protein2ligand-charged',
'bio:active-site:hbond:ligand2protein',
'bio:active-site:hbond:ligand2protein-charged',
'bio:active-site:hbond:non-ideal-protein2ligand',
'bio:active-site:hbond:non-ideal-protein2ligand-charged',
'bio:active-site:hbond:non-ideal-ligand2protein',
'bio:active-site:hbond:non-ideal-ligand2protein-charged',
'bio:active-site:hbond:ligand-intra-molecular',
'bio:active-site:hbond:ligand-intra-molecular-charged',
'bio:active-site:chelator:ligand-chelates',
'bio:active-site:chelator:protein-chelates',
'bio:active-site:chelator:ligand-intra-molecular',
'bio:active-site:chelator:ligand-intra-molecular',
'bio:active-site:salt-bridge:ligand-protein+',
'bio:active-site:salt-bridge:ligand+protein-',
'bio:active-site:stacking:t',
'bio:active-site:stacking:pi',
'bio:active-site:cationpi:ligandpi',
'bio:active-site:cationpi:proteinpi',
'bio:active-site:halogen:ligand-electrophile',
'bio:active-site:halogen:ligand-nucleophile',
'bio:active-site:halogen:protein-electrophile',
'bio:active-site:halogen:protein-nucleophile'
]

