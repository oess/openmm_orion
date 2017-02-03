from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.oedock import FREDDocking
from LigPrepCubes.cubes import SMIRFFParameterization, SetIDTagfromTitle, OEBSinkCube

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SmirffComplexMD")

job.description = """
Starting with SMILE Strings, generate multiconformers using OMEGA,
then dock the multiconformer molecules using FRED,
taking the top scoring pose, parameterize the molecules with a SMIRFF ffxml.
Write them out to an oeb.gz
"""

job.classification = [
    ["OpenEye", "Complex Setup"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ifs", description="File containing SMILES")

omega = OEOmegaConfGen('omega')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')
fred_out = OEBSinkCube('fred_out')
fred_out.set_parameters(suffix='docked')

idtag = SetIDTagfromTitle('idtag')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")
smirff_out = OEBSinkCube('smirff_out')
smirff_out.set_parameters(suffix='smirff')

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein')
complex_out = OEBSinkCube('complex_out')
complex_out.set_parameters(suffix='complex')

md_sim = OpenMMSimulation('md_sim')
md_sim.promote_parameter('complex_mol', promoted_name='complex_mol')
md_sim.promote_parameter('steps', promoted_name='steps')
sim_out = OEBSinkCube('sim_out')
sim_out.set_parameters(suffix='simulation')

cubes = [ifs, omega, fred, idtag, smirff, smirff_out,
        complex_setup, complex_out, md_sim, sim_out]

job.add_cubes(*cubes)

ifs.success.connect(omega.intake)
omega.success.connect(fred.intake)
fred.success.connect(idtag.intake)
idtag.success.connect(smirff.intake)

smirff.success.connect(smirff_out.intake)
smirff.success.connect(complex_setup.intake)

complex_setup.success.connect(complex_out.intake)
complex_setup.success.connect(md_sim.intake)

md_sim.success.connect(sim_out.intake)

if __name__ == "__main__":
    job.run()
