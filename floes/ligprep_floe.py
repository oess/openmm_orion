from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SmirffFragPrep")

job.description = """
Starting with SMILE Strings, generate multiconformers using OMEGA,
then dock the multiconformer molecules using FRED,
taking the top scoring pose, parameterize the molecules with a SMIRFF ffxml.
Write them out to an oeb.gz
"""

job.classification = [
    ["OpenEye", "Ligand Preparation"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ifs", description="File containing SMILES")

omega = OEOmegaConfGen('omega')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")

ofs = OEBSinkCube('ofs')
ofs.set_parameters(suffix='smirff')

job.add_cubes(ifs, omega, fred, idtag, smirff, ofs)
ifs.success.connect(omega.intake)
omega.success.connect(fred.intake)
fred.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)
smirff.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
