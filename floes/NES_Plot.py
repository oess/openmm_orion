from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.FEC.RFEC.cubes import PlotRBFEResults

job = WorkFloe("NES Results Correlation with Exptl DeltaG",
               title="NES Results Correlation with Exptl DeltaG")

job.description = """
NES Results Correlation with Exptl DeltaG
"""

job.classification = [['NES Plot']]
job.uuid = "b1e253a6-0ced-4f67-aee7-0d49b14d16c7"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("RBFE Results Reader", title="RBFE Results Reader")
ifs.promote_parameter("data_in", promoted_name="RBFE_results", title='RBFE Results File',
                      description="RBFE Results File")

plot_res = PlotRBFEResults("RBFEPlot")
plot_res.promote_parameter('lig_exp_file', promoted_name='expt_deltaG', default=None)
plot_res.promote_parameter('symmetrize', promoted_name='symmetrize', default=True)

ofs = DatasetWriterCube('Out', title='Affinities')
ofs.promote_parameter("data_out", promoted_name="out", title="Affinities", description="NES Corr out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="NES Corr Failures out")

job.add_cubes(ifs, plot_res, ofs, fail)

ifs.success.connect(plot_res.intake)
plot_res.success.connect(ofs.intake)
plot_res.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
