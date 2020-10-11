# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
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
import os
import glob

from math import ceil
from uuid import uuid4

from jinja2 import Environment

import pytraj
import mdtraj as md

from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      ComputeCube)

from MDOrion.Standards import Fields

from orionclient.utils import TemporaryPath
from orionclient.session import APISession, in_orion
from orionclient.types import File


from MDOrion.Standards.mdrecord import MDDataRecord


OCCUPANCY_TEMPLATE = """
<html>
<body>
    <div id="ngl-viewer" style="width:100%; height:100%;"></div>
</body>
    <script type="text/javascript" defer="defer">
        loadJavascript = function(url) {
            var xhrObj = new XMLHttpRequest();;
            // open and send a synchronous request
            xhrObj.open('GET', url, false);
            xhrObj.send();
            // Eval the code, to get it into this scope
            // trying to have a <script> block seems to attach it
            // to the parent window, which the iframe doesn't have access
            // to.
            eval(xhrObj.responseText);
        }
        loadJavascript("https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js")

        // Create NGL Stage object
        var stage = new NGL.Stage("ngl-viewer");

        // Handle window resizing
        window.addEventListener( "resize", function( event ){
            stage.handleResize();
        }, false );


        stage.loadFile(
            "https://orion-qa.eyesopen.us/api/v3/storage/files/{{coord_file_id}}/download",
            { ext: "pdb" }
        ).then(function (o) {
          o.addRepresentation("cartoon", { color: "sstruc" })
          o.addRepresentation("licorice", { sele: "not protein and not water"})
          o.autoView()
          {% if traj_file_id %}
          NGL.autoLoad("https://orion-qa.eyesopen.us/api/v3/storage/files/{{traj_file_id}}/download", { ext: "dcd" }).then(function (frames) {
            var traj = o.addTrajectory(frames).trajectory
            var player = new NGL.TrajectoryPlayer(
                traj,
                {
                    step: 1,
                    timeout: 70,
                    start: 0,
                    end: traj.numframes,
                    interpolateType: 'linear',
                    mode: 'once',
                }
            )
              player.play()
          })
          {% endif %}
        });
    </script>
</html>
"""


class OccupancyCalculator(RecordPortsMixin, ParallelMixin, ComputeCube):
    title = "Compute Occupancy"

    def process(self, record, port):
        mdrecord = MDDataRecord(record)

        try:
            title = mdrecord.get_title
        except ValueError:
            title = "Unknown"

        # self.log.info('{}: Attempting MD Traj conversion into OEMols'.format(system_title))

        traj_fn = mdrecord.get_stage_trajectory()

        # self.log.info('{} Temp Directory: {}'.format(system_title, os.path.dirname(traj_fn)))
        # self.log.info('{} Trajectory filename: {}'.format(system_title, traj_fn))

        void, traj_ext = os.path.splitext(traj_fn)

        pdb_fn = None

        if traj_ext == '.h5':
            trj = md.load_hdf5(traj_fn)
        elif traj_ext == '.trr':
            traj_dir = os.path.dirname(traj_fn)
            pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
            trj = md.load_trr(traj_fn, top=pdb_fn)
            trj = trj[1:]
        else:
            raise ValueError("Invalid traj ext: {}".format(traj_ext))
        # Multiple by 10 to go from nm to angstroms
        # bounding = max(*[max(box.x, box.y, box.z) for box in trj.openmm_boxes(0)]) * 10
        # print("BOUNDING", bounding)
        if pdb_fn is None:
            pdb_fn = "random{}".format(uuid4())
            trj[0].save_pdb(pdb_fn)
        coord_file = File.upload(APISession, "PDB Trajectory", pdb_fn)
        APISession.tag_resource(coord_file, "archived")
        # bin_size = 0.5
        # print("BOUNDS", bounding, trj.n_frames)
        # resname = "XIU"  # Figure out if this is always the same
        with TemporaryPath(suffix=".dcd") as temp:
            trj.save_dcd(temp)
            traj_file = File.upload(APISession, "DCD Trajectory", temp)
            APISession.tag_resource(traj_file, "archived")
            # cpp_traj = pytraj.load(temp, pdb_fn)
            # cpp_traj = cpp_traj.superpose(ref=0, mask="@CA")
            # command = "{bounding} {bin_size} {bounding} {bin_size} {bounding} {bin_size} :{resname} normframe".format(
            #     bounding=int(ceil(bounding)),
            #     resname=resname,
            #     bin_size=bin_size
            # )
            # for grid in pytraj.grid(traj=cpp_traj, command=command):
            #     pass
            #     for x, val in enumerate(grid):
            #         for y in range(len(val)):
            #             for z in range(len(val[y])):
            #                 if val[y][z] > 0:
            #                     print("Value", x, y, z, val[y][z])
            #         print(len(val), len(val[0]), len(val[0]), val)
        template = Environment().from_string(OCCUPANCY_TEMPLATE)
        with open("floe_report.html", "w") as ofs:
            ofs.write(template.render(traj_file_id=traj_file.id, coord_file_id=coord_file.id))
        if in_orion():
            report = File.upload(APISession, "{} report".format(title), "floe_report.html")
            tags = ["floe_report"]
            if os.environ["ORION_JOB_ID"]:
                tags.append("Job {}".format(os.environ["ORION_JOB_ID"]))
            APISession.tag_resource(report, *tags)
        self.success.emit(record)
