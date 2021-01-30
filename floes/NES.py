from floe.api import (WorkFloe,
                      ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube


from MDOrion.System.cubes import (CollectionSetting,
                                  ParallelRecordSizeCheck)

from MDOrion.FEC.RFEC.cubes import (BoundUnboundSwitchCube,
                                    RBFECEdgeGathering,
                                    ParallelGMXChimera,
                                    ParallelNESGMX,
                                    NESAnalysis)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import MDFloeReportCube

job = WorkFloe("Non Equilibrium Switching", title="Non Equilibrium Switching")

job.description = """
TBD
"""

job.classification = [['FEC']]
job.uuid = "74cd690f-f98a-47e0-bfa4-1858e4080dc3"
job.tags = [tag for lists in job.classification for tag in lists]

# Unbound and Bound Reader
iunbound_bound = DatasetReaderCube("UnboundBoundReader", title="Unbound amd Bound Reader")
iunbound_bound.promote_parameter("data_in", promoted_name="unbound_bound",
                                 title="Unbound and Bound Input Dataset",
                                 description="Unbound and Bound Dataset")


# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

# Switching Bound and Unbound runs
switch = BoundUnboundSwitchCube("Bound/Unbound Switch", title='Bound/Unbound Switch')

gathering = RBFECEdgeGathering("Gathering", title="Gathering Equilibrium Runs")
gathering.promote_parameter('map_file', promoted_name='map')

chimera = ParallelGMXChimera("GMXChimera", title="GMX Chimera")

unbound_nes = ParallelNESGMX("GMXUnboundNES", title="GMX Unbound NES")
unbound_nes.promote_parameter("time", promoted_name="nes_time", default=0.05)

bound_nes = ParallelNESGMX("GMXBoundNES", title="GMX Bound NES")
bound_nes.promote_parameter("time", promoted_name="nes_time")

# This cube is necessary for the correct working of collections and shards
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

nes_analysis = NESAnalysis("NES_Analysis")

report = MDFloeReportCube("report", title="Floe Report")
report.set_parameters(floe_report_title="NES Report")

check_rec = ParallelRecordSizeCheck("Record Check Success")

ofs = DatasetWriterCube('ofs', title='NES Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="NES Dataset Out",
                      description="NES Dataset Out")

fail = DatasetWriterCube('fail', title='NES Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="NES Failures",
                       description="NES Dataset Failures out")

job.add_cubes(iunbound_bound, coll_open, switch, gathering,
              chimera,  unbound_nes, bound_nes,
              coll_close, nes_analysis, report,
              check_rec, ofs, fail)

# Ligand Setting
iunbound_bound.success.connect(coll_open.intake)
coll_open.success.connect(switch.intake)
switch.success.connect(gathering.intake)
switch.bound_port.connect(gathering.bound_port)

# Chimera NES Setting
gathering.success.connect(chimera.intake)

chimera.success.connect(unbound_nes.intake)
unbound_nes.success.connect(coll_close.intake)

chimera.bound_port.connect(bound_nes.intake)
bound_nes.success.connect(coll_close.intake)

coll_close.success.connect(nes_analysis.intake)
nes_analysis.success.connect(report.intake)
report.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail port connections
coll_open.failure.connect(check_rec.fail_in)
switch.failure.connect(check_rec.fail_in)

gathering.failure.connect(check_rec.fail_in)
chimera.failure.connect(check_rec.fail_in)
unbound_nes.failure.connect(check_rec.fail_in)
bound_nes.failure.connect(check_rec.fail_in)

coll_close.failure.connect(check_rec.fail_in)
nes_analysis.failure.connect(check_rec.fail_in)
report.failure.connect(check_rec.fail_in)
check_rec.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()


