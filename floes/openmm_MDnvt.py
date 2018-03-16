from __future__ import unicode_literals
from floe.api import WorkFloe
from OpenMMCubes.cubes import OpenMMnvtSetCube
from cuberecord import DataSetReaderCube, DataSetWriterCube

from LigPrepCubes.ports import DataSetWriterCubeStripCustom

job = WorkFloe("NVT Simulation")

job.description = """
NVT simulation of an OpenMM-ready System

Ex: python floes/openmm_MDnvt.py --system complex.oeb --ofs-data_out nvt.oeb --picosec 10.0

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature after warming

Outputs:
--------
ofs: Outputs the constant temperature and volume system
"""

job.classification = [['NVT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DataSetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

nvt = OpenMMnvtSetCube('nvt')
nvt.promote_parameter('time', promoted_name='picosec', default=10.0)
nvt.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                      description='Selected temperature in K')
# Restraints
nvt.promote_parameter('restraints', promoted_name='restraints', default='noh (ligand or protein)')
nvt.promote_parameter('restraintWt', promoted_name='restraintWt', default=5.0, description='Restraint weight')

# Trajectory and logging info frequency intervals
nvt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=0.5,
                      description='Trajectory saving interval in ps')
nvt.promote_parameter('reporter_interval', promoted_name='reporter_interval', default=1.0,
                      description='Reporter saving interval in ps')

nvt.promote_parameter('outfname', promoted_name='suffix', default='nvt',
                      description='Equilibration suffix name')
nvt.promote_parameter('tar', promoted_name='tar', default=True)

# ofs = DataSetWriterCubeStripCustom('ofs', title='OFS-Success')
ofs = DataSetWriterCube('ofs', title='OFS-Success')

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, nvt, ofs, fail)
ifs.success.connect(nvt.intake)
nvt.success.connect(ofs.intake)
nvt.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
