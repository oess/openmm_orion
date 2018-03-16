from __future__ import unicode_literals
from floe.api import WorkFloe
from cuberecord import DataSetWriterCube
from ComplexPrepCubes.cubes import SolvationSetCube, ForceFieldSetCube
from LigPrepCubes.cubes import LigandSetChargeCube
from LigPrepCubes.ports import LigandSetReaderCube, DataSetWriterCubeStripCustom
from YankCubes.cubes import YankSolvationFECube
from OpenMMCubes.cubes import OpenMMminimizeSetCube, OpenMMnvtSetCube, OpenMMnptSetCube

job = WorkFloe("Solvation Free Energy")

job.description = """
Solvation Free Energy Calculation of small molecules

Ex. python floes/solvation_free_energy --ligands ligands.oeb
--ofs-data_out fe.oeb

Parameters:
-----------
ligands (file): OEB file of the prepared ligands

Outputs:
--------
ofs: Output file
"""

# *************USER SETTING**************
yank_iteration_per_chunk = 500
chunks = 2
# ***************************************

cube_list = []

job.classification = [['Solvation Free Energy']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = LigandSetReaderCube("Ligands", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")
job.add_cube(iligs)
cube_list.append(iligs)

chargelig = LigandSetChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)
job.add_cube(chargelig)
cube_list.append(chargelig)

solvate = SolvationSetCube("Solvation")
solvate.promote_parameter("density", promoted_name="density", title="Solution density in g/ml", default=1.0,
                          description="Solution Density in g/ml")
solvate.promote_parameter("solvents", promoted_name="solvents", title="Solvent components",
                          default='[H]O[H], ClC(Cl)Cl, CS(=O)C, c1ccccc1',
                          description="Comma separated smiles strings of solvent components")
solvate.promote_parameter("molar_fractions", promoted_name="molar_fractions",
                          title="Molar fractions",
                          default='1.0, 0.0, 0.0, 0.0',
                          description="Comma separated strings of solvent molar fractions")
solvate.promote_parameter('distance_between_atoms', promoted_name='distance_between_atoms', default=2.5)
solvate.promote_parameter("padding_distance", promoted_name="padding_distance", default=11.0,
                          description="The largest dimension (in A) of the solute (along the x, y, or z axis) "
                                      "is determined, and a cubic box of size "
                                      "(largest dimension)+2*padding is used")
job.add_cube(solvate)
cube_list.append(solvate)

ff = ForceFieldSetCube("ForceField")
job.add_cube(ff)
cube_list.append(ff)

# Minimization
minimize = OpenMMminimizeSetCube("Minimize")
minimize.promote_parameter('restraints', promoted_name='m_restraints', default="noh ligand",
                           description='Select mask to apply restarints')
minimize.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                           description='Restraint weight in kcal/(mol A^2')
minimize.promote_parameter('hmr', promoted_name='hmr', default=False,
                           description='Hydrogen Mass Repartitioning')
job.add_cube(minimize)
cube_list.append(minimize)

# NVT Warm-up
warmup = OpenMMnvtSetCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_psec', default=20.0,
                         description='Length of MD run in picoseconds')
warmup.promote_parameter('restraints', promoted_name='w_restraints', default="noh ligand",
                         description='Select mask to apply restarints')
warmup.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0,
                         description='Restraint weight in kcal/(mol A^2')
warmup.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ps')
warmup.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0.0,
                         description='Reporter saving interval in ps')
warmup.promote_parameter('outfname', promoted_name='w_outfname', default='warmup',
                         description='Equilibration suffix name')
warmup.promote_parameter('center', promoted_name='center', default=True)
warmup.promote_parameter('hmr', promoted_name='hmr', default=False,
                         description='Hydrogen Mass Repartitioning')
job.add_cube(warmup)
cube_list.append(warmup)

# NPT Equilibration stage
equil = OpenMMnptSetCube('equil', title='equil')
equil.promote_parameter('time', promoted_name='eq_psec', default=20.0,
                        description='Length of MD run in picoseconds')
equil.promote_parameter('restraints', promoted_name='eq_restraints', default="noh ligand",
                        description='Select mask to apply restraints')
equil.promote_parameter('restraintWt', promoted_name='eq_restraintWt', default=0.1,
                        description='Restraint weight in kcal/(mol A^2')
equil.promote_parameter('trajectory_interval', promoted_name='eq_trajectory_interval', default=0.0,
                        description='Trajectory saving interval in ps')
equil.promote_parameter('reporter_interval', promoted_name='eq_reporter_interval', default=0.0,
                        description='Reporter saving interval in ps')
equil.promote_parameter('outfname', promoted_name='eq_outfname', default='equil',
                        description='Equilibration suffix name')
equil.promote_parameter('hmr', promoted_name='hmr', default=False,
                        description='Hydrogen Mass Repartitioning')
job.add_cube(equil)
cube_list.append(equil)

for i in range(0, chunks):
    solvationfe = YankSolvationFECube("SovationFE"+str(i))
    solvationfe.promote_parameter('iterations', promoted_name='iterations'+str(i),
                                  default=yank_iteration_per_chunk*(i+1))
    solvationfe.promote_parameter('nonbondedCutoff', promoted_name='nonbondedCutoff'+str(i), default=10.0)

    solvationfe.promote_parameter('hmr', promoted_name='hmr'+str(i), default=False,
                                  description='Hydrogen Mass Repartitioning')

    if i == 0:
        solvationfe.promote_parameter('rerun', promoted_name='rerun' + str(i), default=False)
    else:
        solvationfe.promote_parameter('rerun', promoted_name='rerun' + str(i), default=True)

    if i == (chunks - 1):
        solvationfe.promote_parameter('analyze', promoted_name='analyze' + str(i), default=True)

    job.add_cube(solvationfe)
    cube_list.append(solvationfe)

# ofs = DataSetWriterCubeStripCustom('ofs', title='OFS-Success')
ofs = DataSetWriterCube('ofs', title='OFS-Success')
job.add_cube(ofs)
cube_list.append(ofs)

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')
job.add_cube(fail)
cube_list.append(fail)

for i in range(0, len(cube_list)-2):
    print(i)
    cube_list[i].success.connect(cube_list[i + 1].intake)
    if i == len(cube_list) - 3:
        cube_list[i].failure.connect(cube_list[i+2].intake)


if __name__ == "__main__":
    job.run()