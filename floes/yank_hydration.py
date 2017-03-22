from __future__ import unicode_literals
"""
Written by John D. Chodera
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from YankCubes.cubes import YankHydrationCube

"""
Testing via the command line:

python floes/yank_hydration.py --molecules examples/data/freesolv_mini.oeb.gz --success success.sdf --failure failure.sdf --simulation_time 0.001 --nsteps_per_iteration 50

"""

job = WorkFloe("YANK small molecule hydration free energies")

job.description = """
# Compute small molecule hydration free energies using YANK.

This Floe processes the provided molecules, computes free energies of transfer from gas to water, and appends the following SDData properties to the original molecules:

* Hydration free energy in kcal/mol: `DeltaG_yank_hydration`
* Standard error estimate in hydration free energy estimate, in kcal/mol: `dDeltaG_yank_hydration`

*IMPORTANT:* The molecules MUST be charged before feeding them into this workflow.

The following advanced options are available:

* temperature (`float`): Temperature in Kelvin (default: `300.0`)
* pressure (`float`): Pressure in atmospheres (default: `1.0`)
* simulation_time (`float`): Simulation time in ns/replica (default: `0.010`)
* timestep (`float`): Timestep in fs (default: `2.0`)
* yaml_template (`str`): YANK YAML file to use (default template provided)

## About this software

See [http://getyank.org](http://getyank.org) for more information on YANK and the alchemical free energy calculations it supports.

YANK is free (libre) open source software licensed under the [MIT License](https://github.com/choderalab/yank/blob/master/LICENSE).

All source code is available at: [http://github.com/choderalab/yank](http://github.com/choderalab/yank)

YANK is produced by the Chodera lab: [http://choderalab.org](http://choderalab.org).
Please help us make it better by contributing code or funds.

"""

job.classification =[["YANK", "Hydration free energies", "OpenMM", "choderalab"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="molecules", description="Input molecules")

yank_cube = YankHydrationCube('yank_hydration', title = 'Yank for hydration free energies')
for parameter_name in ['solvent', 'temperature', 'pressure', 'nsteps_per_iteration', 'simulation_time', 'timestep', 'verbose']:
    description = yank_cube.parameters()[parameter_name].description
    yank_cube.promote_parameter(parameter_name, promoted_name=parameter_name, description=description)

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
