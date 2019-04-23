import base64

import re

_clus_floe_report_header = """
<html>
<head>
<style>

  .cb-floe-report-wrapper * {
    box-sizing: border-box;
    font-family: 'Helvetica Neue', Helvetica;
  }

  .cb-floe-report__row {
    width:100%;
    display:flex;
  }
  .cb-floe-report__sidebar {
    width: 25%;
  }
  .cb-floe-report__content {
    width: 75%;
  }
  .cb-floe-report__column > * {
    width: 100%;
    height: auto;
  }
  .cb-floe-report__dev-alert {
    background-color: rgb(255, 80, 80);
    color: white;
    padding: .5em;
    font-weight: normal;
  }
  .cb-floe-report-element--analysis {
    /*white-space: pre-wrap;*/
    padding: 0 10px 20px 10px;
    font-size: calc(.5vh + .75vw + .25vmin);
  }

  div.cb-floe-report__analysis-table {
    width: 100%;
  }

  div.cb-floe-report__analysis-table-row {
    display: flex;
    justify-content: space-between;
    padding: 0 5px;
  }

  div.cb-floe-report__analysis-table-row:nth-child(1) {
    font-weight: bold;
    border-bottom: 1px solid black;
  }

  div.cb-floe-report__analysis-table-row > span:nth-child(1) {
    flex-basis: 20%;
    text-align: center;
  }

  div.cb-floe-report__analysis-table-row > span:nth-child(2) {
    flex-basis: 15%;
    text-align: right;
  }

  div.cb-floe-report__analysis-table-row > span:nth-child(3) {
    flex-basis: 50%;
    text-align: left;
    margin-left: auto;
  }

  h2.cb-floe-report-element--header {
    margin-bottom: 0;
    text-align: left;
  }
  h3.cb-floe-report-element--header {
    margin-bottom: 0;
    text-align: center;
  }
  .cb-floe-report__dev-alert small{
    color: white !important;
  }

  /* tab styles */
  div.cb-floe-report__tab-wrapper input {
    display: none;
  }

  div.cb-floe-report__tab-wrapper label {
    display: block;
    float: left;
    padding: 5px 10px;
    cursor: pointer;
    border: 2px solid black;
    margin-right: 2px;
  }

  div.cb-floe-report__tab-wrapper input:checked + label {
    background: black;
    color: white;
    cursor: default;
  }

  div.cb-floe-report__tab-wrapper div.cb-floe-report__tab-content {
    display: none;
    padding: 5px 10px;
    clear: left;
    border: 2px solid black;
  }
"""

_clus_floe_report_header2 = """
</style>

<script type="text/javascript">
  // fix2018aug24 document.addEventListener('DOMContentLoaded', function() {
    var svgs = document.querySelectorAll('.cb-floe-report__tab-content svg');

    // Trigger a click on the buttons in all other SVGs that have a corresponding index
    function syncOtherButtons(svgs, idxSVG, idxBtn) {
      svgs.forEach(function(element, idxA, listA) {
        if (idxA == idxSVG) return false;
        var buttons = element.querySelectorAll('g[id^="tab_button_"]');
        buttons.forEach(function(button, idxB, listB) {
          if (idxB !== idxBtn) return false;
          var evt = new Event('click');
          evt.detail = evt.detail || {};
          evt.detail.isSynthetic = true;
          button.dispatchEvent(evt);
        });
      });
    };

    function listenForClick(e) {
      e.detail = e.detail || {};
      // use a flag to prevent scripted click events from causing infinite recursion!
      if (e.detail && !e.detail.isSynthetic) {
        syncOtherButtons(svgs, idxSVG, idxBtn);
        e.stopPropagation();
      }
    }

    // When a button inside an SVG is clicked, trigger clicks on the others
    svgs.forEach(function(element, idxSVG, listA) {
      var buttons = element.querySelectorAll('g[id^="tab_button_"]');
      buttons.forEach(function(button, idxBtn, listB) {
        button.addEventListener('click', function(e) {
          e.detail = e.detail || {};
          // use a flag to prevent scripted click events from causing infinite recursion!
          if (e.detail && !e.detail.isSynthetic) {
            syncOtherButtons(svgs, idxSVG, idxBtn);
            e.stopPropagation();
          }
        });
      });
    });
  // fix2018aug24 });
</script>
"""

_clus_floe_report_midHtml0 = """

</head>
<body>
<div class="cb-floe-report-wrapper">
  <div class="cb-floe-report__row">
  <div class="cb-floe-report__column cb-floe-report__sidebar">
    {query_depiction}
"""

_clus_floe_report_midHtml1 = """
  </div>

  <div class="cb-floe-report__column cb-floe-report__content">
    <h2> Analysis of Short Trajectory MD </h2>

    <div class="cb-floe-report__tab-wrapper">

"""

_clus_floe_report_midHtml2 = """      </div>
    </div>
  </div>

  <div class="cb-floe-report__row">
    <h2 style="text-align: center; width: 100%"> Ligand Clustering based on Active Site Alignment </h2>
  </div>

  <div class="cb-floe-report__row">
    <div class="cb-floe-report__column cb-floe-report__sidebar">
"""

_clus_floe_report_Trailer = """    </div>

    <div class="cb-floe-report__column cb-floe-report__content">
      <h3 class="cb-floe-report-element--header"> Cluster membership of ligand by Trajectory frame </h3>
      {clusters}
      <h3 class="cb-floe-report-element--header"> RMSD of ligand compared to initial pose, colored by cluster </h3>
      {rmsdInit}
    </div>
  </div>

</div>

</body>
</html>"""

def MakeClusterInfoText(dataDict, rgbVec):
    # Generate text string about Clustering information
    #
    text = []
    nFrames = dataDict['nFrames']
    text.append("""
      <br/>
      <div class="cb-floe-report-element--analysis">""")
    text.append('Clustering by ligand RMSD after alignment by active site C_alphas:\n' )
    text.append('        <br>- Cluster method {}\n'.format( dataDict['ClusterMethod']) )
    text.append('        <br>- Using alpha={:.2f}\n'.format( dataDict['HDBSCAN_alpha']))
    text.append('        <br>- Clustered {} frames\n'.format(nFrames) )
    #
    if dataDict['nClusters']<2:
        text.append('        <br>- Produced {} cluster'.format( dataDict['nClusters']))
    else:
        text.append('        <br>- Produced {} clusters'.format( dataDict['nClusters']))
    nOutliers = dataDict['ClusterVec'].count(-1)
    text.append(' with {:4d} outliers\n'.format( nOutliers))
    text.append('<br><br>Clusters larger than 2%: ' )
    #
    text.append("""
        <div class="cb-floe-report__analysis-table">
          <div class="cb-floe-report__analysis-table-row">
            <span>Cluster</span>
            <span>Size</span>
            <span>Status</span>
          </div>\n
""")

    for i, (count, rgb) in enumerate(zip(dataDict['ClusterCounts'], rgbVec)):
        percent = int(100*count/nFrames)
        if percent < 2:
            continue
        status = 'major'
        if percent < 10:
            status = 'minor (no analysis)'
        text.append("""
          <div class="cb-floe-report__analysis-table-row" style="
                background-color: rgb({r}, {g}, {b} );
                color: white;">
            <span>{clusID}</span>
            <span>{percent}%</span>
            <span>{status}</span>
          </div>\n""".format(clusID=i, percent=percent, status=status,
                             r=rgb[0], g=rgb[1], b=rgb[2]))

    text.append("""
        </div>
      </div>
    """)

    return text


def png_to_data_url(png_data):
    return "<img src='data:image/png;base64," + base64.b64encode(png_data).decode('utf-8') + "'>"


def trim_svg(svg):
    run_init = "\nif (init != null) { try {init() } catch(error) {  true; } };\n"
    svg = re.sub('(<script[^>]*> *<!\[CDATA\[)', '\\1' + run_init, svg)

    idx = svg.find('<svg')
    return svg[idx:]


