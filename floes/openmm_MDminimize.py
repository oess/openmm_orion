from __future__ import unicode_literals
from floe.api import WorkFloe
from OpenMMCubes.cubes import OpenMMminimizeSetCube
from cuberecord import DataSetReaderCube, DataSetWriterCube
from LigPrepCubes.ports import DataSetWriterCubeStripCustom

job = WorkFloe("Minimize")

job.description = """
Minimize an OpenMM-ready solvated complex

Ex: python floes/openmm_prepMDminimize.py --system complex.oeb --ofs-data_out min.oeb --steps 1000`

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
steps (int): Number of MD steps to minimize the system. If 0 until convergence will be reached

Outputs:
--------
ofs: Outputs the minimized system
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DataSetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

min = OpenMMminimizeSetCube('Minimize')
min.promote_parameter('steps', promoted_name='steps')

# fs = DataSetWriterCubeStripCustom('ofs', title='OFS-Success')
ofs = DataSetWriterCube('ofs', title='OFS-Success')

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, min, ofs, fail)
ifs.success.connect(min.intake)
min.success.connect(ofs.intake)
min.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
