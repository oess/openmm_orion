from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from LigPrepCubes.omega import OEOmegaConfGen
from LigPrepCubes.cubes import Attachffxml, FREDDocking
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("FRED")

job.description = """
Generate multiple conformers from SMILE Strings and dock with FRED

"""

job.classification = [
    ["OpenEye", "Ligand Preparation"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ifs", description="File containing SMILES")

ffxml = Attachffxml('ffxml')
ffxml.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")

omega = OEOmegaConfGen('omega')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

#complex_setup = OpenMMComplexSetup("complex_setup")
#complex_setup.promote_parameter('protein', promoted_name='protein')
#complex_setup.promote_parameter('molecule_forcefield', promoted_name='molecule_forcefield')
#complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_forcefield')
#complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_forcefield')

#md_sim = OpenMMSimulation('md_sim')
#md_sim.promote_parameter('complex_mol', promoted_name='complex_mol')
#md_sim.set_parameters(complex_mol='complex.oeb.gz')

ofs = OEMolOStreamCube('ofs')
ofs.set_parameters(data_out="test.oeb.gz")


job.add_cubes(ifs, ffxml, omega, fred)
ifs.success.connect(ffxml.intake)
ffxml.success.connect(omega.intake)
omega.success.connect(fred.intake)
#fred.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
