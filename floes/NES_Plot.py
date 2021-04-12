from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.FEC.RFEC.cubes import PlotRBFEResults

job = WorkFloe("Testing NES Plot Results",
               title="NES Results")

job.description = """
Testing NES Plot Results
"""

job.classification = [['NES Plot']]
job.uuid = "b1e253a6-0ced-4f67-aee7-0d49b14d16c7"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='RBFE Input File',
                      description="System input file")

plot_res = PlotRBFEResults("RBFEPlot")
plot_res.promote_parameter('lig_exp_file', promoted_name='exp', default=None)
plot_res.promote_parameter('symmetrize', promoted_name='symmetrize', default=False)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(ifs, plot_res, fail)

ifs.success.connect(plot_res.intake)
plot_res.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
