from __future__ import unicode_literals
from floe.api import WorkFloe

from cuberecord import DataSetWriterCube, DataSetReaderCube
from LigPrepCubes.ports import LigandSetReaderCube
from LigPrepCubes.cubes import LigandSetChargeCube
from LigPrepCubes.ports import DataSetWriterCubeStripCustom

from ComplexPrepCubes.port import ProteinSetReaderCube

from ComplexPrepCubes.cubes import (
    ComplexSetPrepCube,
    HydrationSetCube,
    SolvationSetCube,
    ForceFieldSetCube)

job = WorkFloe("Complex Preparation")

job.description = """
Complex Preparation Workflow

Ex. python floes/openmm_complex_prep.py --protein protein.oeb
--ligands ligands.oeb  --ofs-data_out complex.oeb

Parameters:
-----------
protein (file): OEB file of the prepared protein
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = LigandSetReaderCube("Ligand Reader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigandSetChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

iprot = ProteinSetReaderCube("Protein Reader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title="Protein Input File", description="Protein file name")
iprot.promote_parameter("protein_prefix", promoted_name="protein_prefix", default='PRT',
                        description="Protein Prefix")

complx = ComplexSetPrepCube("Complex")

solvate = HydrationSetCube("Hydration")

ff = ForceFieldSetCube("ForceField")

ofs = DataSetWriterCube('ofs', title='OFS-Success')
#ofs = DataSetWriterCubeStripCustom('ofs', title='OFS-Success')


fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iligs, chargelig, iprot, complx, solvate, ff, ofs, fail)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
iprot.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(ofs.intake)
ff.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
