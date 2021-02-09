from floe.api import WorkFloe

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


# Unbound Reader
iun = DatasetReaderCube("UnboundReader", title="Unbound Reader")
iun.promote_parameter("data_in", promoted_name="unbound", title="Unbound Input Dataset", description="Unbound Input Dataset")

ibn = DatasetReaderCube("BoundReader", title="Bound Reader")
ibn.promote_parameter("data_in", promoted_name="bound", title="Bound Input Dataset", description="Bound Input Dataset")

# Switching Bound and Unbound runs
switch = BoundUnboundSwitchCube("Bound/Unbound Switch", title='Bound/Unbound Switch')

# This cube is necessary for the correct work of collection and shard
coll_open_write = CollectionSetting("OpenWriteNESCollection", title="OpenWriteNESCollection")
coll_open_write.set_parameters(open=True)
coll_open_write.set_parameters(write_new_collection='NES_OPLMD')

gathering = RBFECEdgeGathering("Gathering", title="Gathering Equilibrium Runs")
gathering.promote_parameter('map_file', promoted_name='map')

chimera = ParallelGMXChimera("GMXChimera", title="GMX Chimera")
chimera.promote_parameter("trajectory_frames", promoted_name="trajectory_frames", default=80,
                          description="The total number of trajectory frames to be used along the NE switching")

unbound_nes = ParallelNESGMX("GMXUnboundNES", title="GMX Unbound NES")
unbound_nes.promote_parameter("time", promoted_name="nes_time", default=0.05)

bound_nes = ParallelNESGMX("GMXBoundNES", title="GMX Bound NES")
bound_nes.promote_parameter("time", promoted_name="nes_time")

nes_analysis = NESAnalysis("NES_Analysis")

# This cube is necessary for the correct working of collections and shards
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

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

job.add_cubes(iun, ibn, coll_open_write, switch, gathering,
              chimera,  unbound_nes, bound_nes,
              nes_analysis, coll_close, report,
              check_rec, ofs, fail)

# Readers Setting
iun.success.connect(coll_open_write.intake)
ibn.success.connect(coll_open_write.intake)

coll_open_write.success.connect(switch.intake)

switch.success.connect(gathering.intake)
switch.bound_port.connect(gathering.bound_port)

# Chimera NES Setting
gathering.success.connect(chimera.intake)

chimera.success.connect(unbound_nes.intake)
unbound_nes.success.connect(nes_analysis.intake)

chimera.bound_port.connect(bound_nes.intake)
bound_nes.success.connect(nes_analysis.intake)

nes_analysis.success.connect(report.intake)
report.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail port connections
coll_open_write.failure.connect(check_rec.fail_in)
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


