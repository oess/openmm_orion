from __future__ import unicode_literals
"""
Written by John D. Chodera
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from YankCubes.cubes import YankBindingCube

"""
Testing via the command line:

python floes/yank_binding.py --receptor examples/data/T4-protein.pdb --molecules examples/data/p-xylene.mol2 --success success.sdf --failure failure.sdf --simulation_time 0.0001 --nsteps_per_iteration 5 --minimize 0

"""

job = WorkFloe("YANK small molecule absolute binding free energies")

job.description = """
# Compute small molecule absolute binding free energies using YANK.

This Floe processes the provided molecules, computes free energies of binding to a specified receptor, and appends the following SDData properties to the original molecules:

* Absolute binding free energy in kcal/mol: `DeltaG_yank_binding`
* Standard error estimate in binding free energy estimate, in kcal/mol: `dDeltaG_yank_binding`

*IMPORTANT:* The molecules MUST be charged before feeding them into this workflow.

The following advanced options are available:

* temperature (`float`): Temperature in Kelvin (default: `300.0`)
* pressure (`float`): Pressure in atmospheres (default: `1.0`)
* simulation_time (`float`): Simulation time in ns/replica (default: `0.010`)
* timestep (`float`): Timestep in fs (default: `2.0`)
* yaml_template (`str`): YANK YAML file to use (default template provided)
* minimize (`bool`): If `True`, the complex will be minimized prior to simuation (recommended)

## WARNINGS

This floe is **highly experimental** and currently uses harmonic restraints to ensure the ligand remains close to the receptor.
In principle, this procedure can dock the ligand into the receptor, but in practice, this may not work well for ligands with multiple rotatable bonds.
This has only been tested on small fragment-like ligands so far.

**Use at your own risk!**

## About this software

See [http://getyank.org](http://getyank.org) for more information on YANK and the alchemical free energy calculations it supports.

YANK is free (libre) open source software licensed under the [MIT License](https://github.com/choderalab/yank/blob/master/LICENSE).

All source code is available at: [http://github.com/choderalab/yank](http://github.com/choderalab/yank)

YANK is produced by the Chodera lab: [http://choderalab.org](http://choderalab.org).
Please help us make it better by contributing code or funds.

"""

job.classification =[["YANK", "Binding free energies", "OpenMM", "choderalab"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="molecules", description="Input molecules")

yank_cube = YankBindingCube('yank_binding', title = 'Yank for binding free energies')
for parameter_name in ['receptor', 'solvent', 'temperature', 'pressure', 'nsteps_per_iteration', 'simulation_time', 'timestep', 'minimize', 'verbose']:
    promoted_name = parameter_name
    description = yank_cube.parameters()[parameter_name].description
    yank_cube.promote_parameter(parameter_name, promoted_name=promoted_name, description=description)

success_ofs = OEMolOStreamCube("success_ofs")
success_ofs.promote_parameter("data_out", promoted_name="success", description="Output molecules")

failure_ofs = OEMolOStreamCube("failure_ofs")
failure_ofs.promote_parameter("data_out", promoted_name="failure", description="Failed molecules")

cubes = [ifs, yank_cube, success_ofs, failure_ofs]

job.add_cubes(*cubes)

ifs.success.connect(yank_cube.intake)
yank_cube.success.connect(success_ofs.intake)
yank_cube.failure.connect(failure_ofs.intake)


if __name__ == "__main__":
    job.run()
