#!/bin/env python

from openeye import oechem

ifs = oechem.oemolistream('freesolv_sybyl.oeb.gz')
ofs = oechem.oemolostream('freesolv_mini.oeb.gz')

molecules = [
    'octafluorocyclobutane',
    '2,2-dimethylpentane',
    'but-1-ene',
    'toluene',
    '1-butoxybutane',
    'butan-1-ol',
    '4-methylpyridine',
    'methanol',
    '3-nitroaniline',
    'N-methylacetamide',
    '4-nitrophenol',
    '5-iodouracil'
    ]

mol = oechem.OEGraphMol()
for mol in ifs.GetOEGraphMols():
    if mol.GetTitle() in molecules:
        print(mol.GetTitle())
        oechem.OEWriteMolecule(ofs, mol)

ofs.close()
ifs.close()
