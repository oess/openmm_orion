# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from os import path

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube


from MDOrion.System.cubes import (ParallelRecordSizeCheck)

from floes.SubfloeFunctions import (setup_NonEquilSwch_GMX)


job = WorkFloe("Non-Equilibrium Switching", title="Non-Equilibrium Switching")

job.description = open(path.join(path.dirname(__file__), 'NonEquilSwch.rst'), 'r').read()

job.classification = [['FEC']]
job.uuid = "74cd690f-f98a-47e0-bfa4-1858e4080dc3"
job.tags = [tag for lists in job.classification for tag in lists]


# Unbound Reader
iun = DatasetReaderCube("UnboundReader", title="Unbound Reader")
iun.promote_parameter("data_in", promoted_name="unbound", title="Unbound Input Dataset", description="Unbound Input Dataset")

ibn = DatasetReaderCube("BoundReader", title="Bound Reader")
ibn.promote_parameter("data_in", promoted_name="bound", title="Bound Input Dataset", description="Bound Input Dataset")

check_rec = ParallelRecordSizeCheck("Record Check Success")

job.add_cubes(iun, ibn, check_rec)

# Call subfloe function to start up the MD, equilibrate, and do the production run
NonEquilSwch_startup_options = {}
NonEquilSwch_startup_options['edge_map_file'] = 'map'
NonEquilSwch_startup_options['n_traj_frames'] = 80
NonEquilSwch_startup_options['nes_switch_time_in_ns'] = 0.05
setup_NonEquilSwch_GMX(job, ibn, iun, check_rec, NonEquilSwch_startup_options)


if __name__ == "__main__":
    job.run()


