import os
import traceback
import re


from floe.api import (ParallelMixin,
                      parameter)

# Just for old orion testing
from datarecord import (Types,
                        Meta,
                        OEFieldMeta,
                        OEField,
                        OERecord)


from cuberecord import OERecordComputeCube
from Standards import (Fields,
                       MDRecords,
                       MDStageNames)

from orionclient.session import in_orion, OrionSession
from orionclient.types import File
from os import environ

from openeye import oechem, oedepict


_clus_floe_report_template = """
<html>
<head>
<style>

  .row {{
    width:100%;
    display:flex;
  }}
  .sidebar {{
    width: 25%;
  }}
  .content {{
    width: 75%;
  }}
  .column > * {{
    width: 100%;
    height: auto;
  }}
  pre {{
    white-space: pre-wrap;
    font-family: 'Helvetica Neue', Helvetica;
    padding: 0 10px 20px 10px;
  }}
  h2 {{
    margin-bottom: 0;
    text-align: left;
  }}
  h3 {{
    margin-bottom: 0;
    text-align: center;
  }}
</style>

</head>
<body>

<div class="row">
<div class="column sidebar">
  {query_depiction}
  {rmsd_hist}
</div>

<div class="column content">
  <h2> Analysis of Short Trajectory MD </h2>
  {traj}
</div>
</div>

<div class="row">
  <h2 style="text-align: center; width: 100%"> Ligand Clustering based on Active Site Alignment </h2>
</div>

<div class="row">
<div class="column sidebar">
  <br/><br/>
  <pre>
  {analysis}
  </pre>
</div>

<div class="column content">
  <h3> Cluster membership of Trajectory ligand by Trajectory frame </h3>
  {clusters}
  <h3> RMSD of ligand compared to initial pose, colored by cluster </h3>
  {rmsdInit}
</div>
</div>

</body>
</html>
"""

def _trim_svg(svg):
    # svg = svg.decode('utf-8')
    # insert a call to init.. sigh
    run_init = "\nif (init != null) { try {init() } catch(error) {  true; } };\n" 
    svg = re.sub('(<script[^>]*> *<!\[CDATA\[)', '\\1' + run_init, svg)

    idx = svg.find('<svg')
    return svg[idx:]

def CheckAndGetValueFull( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))

def CheckAndGetValue( record, field):
    if not record.has_value(field):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field.get_name() ))
        raise ValueError('The record does not have field {}'.format( field.get_name() ))
    else:
        return record.get_value(field)


#class MDTrajAnalysisClusterReport(ParallelMixin, OERecordComputeCube):
class MDTrajAnalysisClusterReport(OERecordComputeCube):
    title = 'Extract relevant outputs of MD Traj Cluster  Analysis'

    version = "0.1.0"
    classification = [["Simulation", "Traj Analysis"]]
    tags = ['Parallel Cube']

    description = """
    Extract relevant outputs of Ligand and Protein
    Short Traj MD Traj Analysis and write them to files.

    This cube takes as input the OERecord containing the work
    product of trajectory analysis on Short Traj MD results.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # title of entire solvated protein-ligand system
            opt['Logger'].info('Starting Floe Report generation for MD Traj Analysis')
            system_title = CheckAndGetValue( record, Fields.title)
            opt['Logger'].info('{} Attempting to extract MD Traj Analysis results'
                .format(system_title) )
            ligInitPose = CheckAndGetValue( record, Fields.ligand)

            # Extract the traj SVG from the OETraj record
            analysesDone = CheckAndGetValueFull( record, 'AnalysesDone', Types.StringVec)
            if 'OETraj' not in analysesDone:
                raise ValueError('{} does not have OETraj analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found OETraj analyses'.format(system_title) )
            # Extract the relevant traj SVG from the OETraj record
            oetrajRecord = CheckAndGetValueFull( record, 'OETraj', Types.Record)
            opt['Logger'].info('{} found OETraj record'.format(system_title) )
            trajSVG = CheckAndGetValueFull( oetrajRecord, 'TrajSVG', Types.String)

            # Extract the three plots from the TrajClus record
            analysesDone = CheckAndGetValueFull( record, 'AnalysesDone', Types.StringVec)
            if 'TrajClus' not in analysesDone:
                raise ValueError('{} does not have TrajClus analyses done'.format(system_title) )
            else:
                opt['Logger'].info('{} found TrajClus analyses'.format(system_title) )
            # Extract the relevant traj SVG from the TrajClus record
            clusRecord = CheckAndGetValueFull( record, 'TrajClus', Types.Record)
            opt['Logger'].info('{} found TrajClus record'.format(system_title) )
            trajHistRMSD_svg = CheckAndGetValueFull( clusRecord, 'HistSVG', Types.String)
            trajClus_svg = CheckAndGetValueFull( clusRecord, 'ClusSVG', Types.String)
            rmsdInit_svg = CheckAndGetValueFull( clusRecord, 'rmsdInitPose', Types.String)
            #trajHeatRMSD_png = CheckAndGetValueFull( clusRecord, 'HeatPNG', Types.Blob)
            opt['Logger'].info('{} found the TrajClus plots'.format(system_title) )

            # Generate text string about Clustering information
            analysis_txt = []
            nFrames = CheckAndGetValueFull(clusRecord, 'nFrames', Types.Int)
            analysis_txt.append('Clustering {} frames\n'.format( nFrames))
            clusMethod = CheckAndGetValueFull(clusRecord, 'ClusterMethod', Types.String)
            alpha = CheckAndGetValueFull(clusRecord, 'HDBSCAN_alpha', Types.Float)
            analysis_txt.append('Cluster method {} with alpha {:.2f}\n'.format( clusMethod, alpha))
            nClusters = CheckAndGetValueFull(clusRecord, 'nClusters', Types.Int)
            analysis_txt.append('produced {} clusters:\n'.format( nClusters))
            clusCounts = CheckAndGetValueFull(clusRecord, 'ClusterCounts', Types.IntVec)
            for i, count in enumerate(clusCounts):
                analysis_txt.append('cluster {} contains {} frames\n'.format( i, count))

            opt['Logger'].info('{} finished writing analysis files'.format(system_title) )

            oedepict.OEPrepareDepiction(ligInitPose)
            img = oedepict.OEImage(400, 300)
            oedepict.OERenderMolecule(img, ligInitPose)

            with open("md_clus_report.html", 'a') as report_file:
                report_file.write(_clus_floe_report_template.format(
                    query_depiction=oedepict.OEWriteImageToString("svg", img).decode("utf8"),
                    rmsd_hist=_trim_svg(trajHistRMSD_svg),
                    analysis="".join(analysis_txt),
                    clusters=_trim_svg(trajClus_svg),
                    traj=_trim_svg(trajSVG),
                    rmsdInit=_trim_svg(rmsdInit_svg),
                    )
                )

                report_file.close()

            if in_orion():
                session = OrionSession()

                file_upload = File.upload(session, "{} MD Cluster Report".format(system_title), "./md_clus_report.html")
                session.tag_resource(file_upload, "floe_report")
                job_id = environ.get('ORION_JOB_ID')
                if job_id:
                    session.tag_resource(file_upload, "Job {}".format(job_id))

            self.success.emit(record)

        except:
            print(traceback.format_exc())
            #self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

