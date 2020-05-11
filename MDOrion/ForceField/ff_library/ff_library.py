import  os

import glob

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
COFACTOR_DIR = os.path.join(PACKAGE_DIR, "MDOrion/ForceField/ff_library", "cofactors")

proteinff = {'Amber99SBildn': 'amber99sbildn.xml',
             'Amber99SB': 'amber99sb.xml',
             'Amber14SB': 'amber14/protein.ff14SB.xml',
             'AmberFB15': 'amberfb15.xml'}

protein_extended_ff = {'Amber14SB': os.path.join(PACKAGE_DIR, 'MDOrion/ForceField/ffext/amber14SB_extendend.xml')}

ligandff = {'Gaff': 'GAFF',
            'Gaff2': 'GAFF2',
            'Smirnoff99Frosst': 'smirnoff99Frosst.offxml',
            'OpenFF_1.0.0': "openff_unconstrained-1.0.0.offxml",
            'OpenFF_1.1.0': "openff_unconstrained-1.1.0.offxml"}

otherff = {'Gaff': 'GAFF',
           'Gaff2': 'GAFF2',
           'Smirnoff99Frosst': 'smirnoff99Frosst.offxml',
           'OpenFF_1.0.0': "openff_unconstrained-1.0.0.offxml",
           'OpenFF_1.1.0': "openff_unconstrained-1.1.0.offxml"}


solventff = {'Tip3p': 'tip3p.xml'}

counter_ionsff = {'Counter_ions': 'amber14/tip3p.xml'}

metals_ff = {'Metals': 'amber14/tip3p.xml'}

excipients_ff = {'Excipients': 'amber14/tip3p.xml'}

cofactors_ff = {'Cofactors': [f for f in glob.glob(COFACTOR_DIR+"/*.xml")]}

lipids_ff = {'Lipids': 'amber14/lipid17.xml'}

nucleics_ff = {'Nucleics': 'amber14-all.xml'}


protein_standard_residue_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                  'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                  'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                  'SER', 'THR', 'TRP', 'TYR', 'VAL',
                                  'ASX', 'GLX', 'CYX', 'CYH', 'HID',
                                  'HIE', 'HIP']
