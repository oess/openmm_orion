from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization, FREDDocking

job = WorkFloe("SmilesLigPrep")

job.description = """
**Prepare a molecule dataset from SMILES strings**

Parse SMILES, generate multiconformer OEMolecules, assign partial charges, dock using
the FRED docking engine and parameterize the ligand with the supported forcefield parameters
(GAFF/GAFF2/SMIRNOFF).
"""

job.classification = [["Ligand Preparation"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="File containing SMILES")

charge = ChargeMCMol('charge')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

lig_param = LigandParameterization('lig_param')
lig_param.promote_parameter('molecule_forcefield', promoted_name='molecule_forcefield', description='Forcefield for molecule')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')

job.add_cubes(ifs, charge, fred, lig_param, ofs, fail)
ifs.success.connect(charge.intake)
charge.success.connect(fred.intake)
fred.success.connect(lig_param.intake)
lig_param.success.connect(ofs.intake)
lig_param.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
