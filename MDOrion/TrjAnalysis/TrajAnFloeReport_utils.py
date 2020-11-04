import base64

import re

import oetrajanalysis.Clustering_utils as clusutl

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
  .cb-floe-report__plot {
    padding-left: 4vw;
    padding-right: 4vw;
    text-align: center;
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
    text-align: center;
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
  <br><hr>
"""

_clus_floe_report_midHtml2a = """      </div>
    </div>
  </div>
  <br><hr>
  <div class="cb-floe-report__row">
    <h2 style="text-align: center; width: 100%"> Ligand Clustering based on Active Site Alignment </h2>
    <br>
  </div>

  <div class="cb-floe-report__row">
"""
#  <div class="cb-floe-report__row">
#    <div class="cb-floe-report__column cb-floe-report__sidebar">
#"""

_clus_floe_report_stripPlots = """    </div>

    <div class="cb-floe-report__column cb-floe-report__plot">
      <h2 class="cb-floe-report-element--header"> Cluster membership of ligand by Trajectory frame </h2>
      {clusters}
    </div>"""
#    </div>
#  </div>"""
# removed from _clus_floe_report_stripPlots but commented out here in case needed again
#      <h3 class="cb-floe-report-element--header"> RMSD of ligand compared to initial pose, colored by cluster </h3>
#      {rmsdInit}


_clus_floe_report_Trailer = '\n</body>\n</html>\n'

def MakeClusterInfoText(dataDict, popResults, rgbVec):
    # Generate text string about Clustering information
    #
    text = []
    nFrames = dataDict['nFrames']
    nMajor = dataDict['nMajorClusters']
    nMajorThresh = dataDict['MajorClusThreshold']
    rowColors = clusutl.ColorblindHexMarkerColors(nMajor) + ['#4c4c4c']
    #
    text.append("""
      <br/>
      <div class="cb-floe-report-element--analysis">""")
    text.append('Clustering by ligand RMSD after alignment by active site C_alphas:\n' )
    text.append('        <br>- Cluster method {}\n'.format( dataDict['ClusterMethod']) )
    #text.append('        <br>- Using alpha={:.2f}\n'.format( dataDict['HDBSCAN_alpha']))
    text.append('        <br>- Clustered {} frames\n'.format(nFrames) )
    #
    if dataDict['nClusters']<2:
        text.append('        <br>- Produced {} cluster'.format( dataDict['nClusters']))
    else:
        text.append('        <br>- Produced {} clusters'.format( dataDict['nClusters']))
    nOutliers = dataDict['ClusterVec'].count(-1)
    text.append(' and {:4d} Outliers,'.format( nOutliers))
    text.append(' with {} Major clusters (>{:.0%} of trajectory)\n'.
                format( nMajor,nMajorThresh) )

    text.append('<br><br>Major Clusters (> {:.0%} of trajectory):<br>'.format( nMajorThresh) )
    #
    text.append("""
        <div class="cb-floe-report__analysis-table">
          <div class="cb-floe-report__analysis-table-row">
            <span>Cluster</span>
            <span>Size</span>
            <span>MMPBSA* &plusmn; StdErr</span>
          </div>\n""")
           #<span>&lt;MMPBSA&gt;* &plusmn; StdErr</span>

    clusSize = popResults['ClusTot']
    ByClusMean = popResults['OEZap_MMPBSA6_ByClusMean']
    ByClusSerr = popResults['OEZap_MMPBSA6_ByClusSerr']
    for i, (count, color) in enumerate(zip(clusSize, rowColors)):
        percent = count/nFrames
        mmpbsa = '{:6.2f} &plusmn; {:4.2f}'.format(ByClusMean[i], 1.96 * ByClusSerr[i])
        if i<nMajor:
            clusName = str(i)
        else:
            clusName = 'Other**'
        text.append("""
          <div class="cb-floe-report__analysis-table-row" style="
                background-color: {hexColor};
                color: white;">
            <span>{clusID}</span>
            <span>{percent:.1%}</span>
            <span>{mmpbsadata}</span>
          </div>\n""".format(clusID=clusName, percent=percent, mmpbsadata=mmpbsa,
                             hexColor=color))
    # Footnotes to Table.
    text.append('*Values in kcal/mol.\n')
    text.append('**Other includes Minor clusters as well as Outliers from the clustering.\n')

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


def HtmlGridColColorRule(rulePrefix,nColums,colToColor,color):
    stride = ':nth-child({:d}n+{:d})'.format(nColums,colToColor)
    rule = rulePrefix+stride+' {\n  background-color: '+color+';\n  color: white;\n}\n\n'
    return rule


def HtmlGridRowColorRule(rulePrefix,nColums,rowToColor,color):
    strt = 1+(rowToColor-1)*nColums
    end = rowToColor*nColums
    stride = ':nth-child(n+{:d}):nth-child(-n+{:d})'.format(strt, end)
    rule = rulePrefix+stride+' {\n  background-color: '+color+';\n  color: white;\n}\n\n'
    return rule


def HtmlTableType1(tabName,colNames,rowNames,tabVals,colors,colorBy='column',fieldWidth=12):
    # Setup
    nCols = len(colNames)
    nRows = len(rowNames)+1
    #
    # Header
    grdHeader = tabName+'-container {\n  display: grid;\n  grid-template-columns:'
    for i in range(nCols):
        grdHeader += ' auto'
    grdHeader += ';\n  grid-template-rows:'
    for i in range(nRows):
        grdHeader += ' 1fr'
    grdHeader += ';\n  grid-row-gap: 10px;\n  grid-column-gap: 0px;\n  padding: 2px;\n'
    grdHeader += '  max-width: '+str(15*fieldWidth*nCols)+'px;\n'
    grdHeader += '  margin: auto;\n'
    grdHeader += '  background-color: #a0a0a0;\n}\n'
    grdHeader += '\n.'+tabName+'-bold {\n  font-weight: bold;\n}\n'
    #
    # General rule for elements in the table
    grdItemRule = '\n'+tabName+'-item  {\n  background-color: #fafafa;\n'
    grdItemRule += '  text-align: center;\n  padding: 10px 0;\n  font-size: 20px;\n}\n'
    #
    # Specific child rules to color columns or rows according to color vector
    grdChildRules = ''
    for i,color in enumerate(colors):
        if color is None:
            continue
        if colorBy=='column':
            grdChildRules += HtmlGridColColorRule(tabName+'-item',nCols,i+1,color)
        else:
            grdChildRules += HtmlGridRowColorRule(tabName+'-item',nCols,i+2,color)
    #
    grdStyle = grdHeader+grdItemRule+grdChildRules
    #
    # table body
    strtBold = '  <'+tabName+'-item class="STMD_ClusFracByConf-bold">'
    strt = '  <'+tabName+'-item>'
    end = '</'+tabName+'-item>\n'
    grdBody = '\n<'+tabName+'-container>\n'
    for field in colNames:
        grdBody += strtBold+field+end
    for name, rowVals in zip(rowNames, tabVals):
        grdBody += strtBold+name+end
        for val in rowVals:
            grdBody += strt+val+end
    grdBody+= '</'+tabName+'-container>\n\n'
    #
    return grdStyle, grdBody


def HtmlWriteTableBody(htmlTableDict,headermods=''):
    html = ''
    chunk = htmlTableDict['title']
    if chunk is not None:
        html += '  <br><br>\n  <h2 class="'+headermods+'" style="text-align: center;">'+chunk+'</h2><br>\n'
    chunk = htmlTableDict['subtitle']
    if chunk is not None:
        html += '  <h3 class="'+headermods+'">'+chunk+'</h3><br>\n'
    chunk = htmlTableDict['body']
    if chunk is not None:
        html += chunk
    chunk = htmlTableDict['descr']
    if chunk is not None:
        html += '  <br>\n  <p class="'+headermods+'">'+chunk+'</p><br>\n'
    html += '<hr style="max-width: 800px;">\n'
    return html


def HtmlMakeClusterPopTables(popResults):
    '''Generate the html for the Population analysis of the Major Chlusters, the results
    of which are passed in as dictionaries (the three input arguments).
    Three tables will be generated: the first for the general information on Clusters,
    the second for the cluster occupancy of each starting pose's MD trajectory,
    and the third for the RMSD deviation of each major cluster from each starting pose.'''
    grdPrefix = 'STMD_'
    nConfs = popResults['nConfs']
    nMajorPlus1 = popResults['nMajorPlus1']

    # Column colors are None for col 1, then cluster colors to last col, then gray for Other
    colorBy = 'column'
    rowColors = clusutl.ColorblindHexMarkerColors(nMajorPlus1 - 1) + ['#4c4c4c']
    colColors = [None] + rowColors
    # print(colColors)

    # Table for Cluster size % of overall sampling and ensemble MMPBSA values
    tabName = 'ClusInfo'
    grdName = grdPrefix + tabName
    ByClusMean = popResults['OEZap_MMPBSA6_ByClusMean']
    ByClusSerr = popResults['OEZap_MMPBSA6_ByClusSerr']
    clusPop = popResults['ClusTot']
    # calculate cluster % occupancy in whole trajectory
    clusPercent = [clus/popResults['MaxIdx'] for clus in clusPop]
    # calculate MMPBSA mean and stderr for each cluster
    rows = len(clusPop)
    strVals = []
    for row in range(rows):
        vals = [str(clusPop[row]), '{:.1%}'.format(clusPercent[row])]
        vals.append('{:6.2f} &plusmn; {:4.2f} kcal/mol'.format(ByClusMean[row], 1.96 * ByClusSerr[row]))
        strVals.append(vals)
    #print(strVals)
    # format column and row names
    colNames = ['Cluster', 'Size', '% of Trajectory', '&lt;MMPBSA&gt; &plusmn; StdErr']
    rowNames = ['Cluster ' + str(i) + ' :' for i in range(len(clusPop))]
    rowNames[-1] = "Other :"
    # assemble dict of html components to go into the table formatting function
    html_ClusInfo = dict()
    html_ClusInfo['title'] = 'Cluster Information'
    html_ClusInfo['subtitle'] = None
    colorBy = 'row'
    maxFieldChars = 10
    style, body = HtmlTableType1(grdName, colNames, rowNames, strVals, rowColors, colorBy, maxFieldChars)
    html_ClusInfo['style'] = style
    html_ClusInfo['body'] = body
    descr = 'This Table gives basic information on each major cluster from the whole trajectory: '
    descr += 'how large is the cluster in terms of number of frames, percent of the overall '
    descr += 'trajectory, and the ensemble average MMPBSA interaction energy (with standard error). '
    descr += 'The "Other" category includes minor clusters as well as Outliers from the clustering.'
    html_ClusInfo['descr'] = descr

    # Table for Cluster population % by traj for each starting Pose
    tabName = 'ClusFracByConf'
    grdName = grdPrefix + tabName
    vals = popResults['ClusFracByConf']
    # print(vals)
    strVals = [['{:.0%}'.format(vals[i][j]) for j in range(nMajorPlus1)] for i in range(nConfs)]
    # format column and row names
    nPoseToPoseID = popResults['nPoseToPoseID']
    poseNames = ['Pose ' + str(nPoseToPoseID[i]) + ' :' for i in range(nConfs)]
    # print(poseNames)
    colNames = ['Pose']
    for i in range(nMajorPlus1):
        colNames.append('Cluster ' + str(i))
    colNames[-1] = "Other"
    # assemble dict of html components to go into the table formatting function
    html_ClusFracByConf = dict()
    html_ClusFracByConf['title'] = 'Cluster Percentage by Starting Pose'
    html_ClusFracByConf['subtitle'] = "How much time did each starting pose's trajectory spend in each major cluster?"
    maxFieldChars = 9
    colorBy = 'column'
    style, body = HtmlTableType1(grdName, colNames, poseNames, strVals, colColors, colorBy, maxFieldChars)
    html_ClusFracByConf['style'] = style
    html_ClusFracByConf['body'] = body
    descr = 'For each starting pose for an MD run, this Table shows the proportion '
    descr += 'of time the ligand spent in each major cluster. '
    descr += 'The "Other" category includes all minor clusters as well as Outliers from the clustering.'
    html_ClusFracByConf['descr'] = descr

    # Table for Cluster RMSDs to starting Poses
    tabName = 'ConfRMSDsByClus'
    grdName = grdPrefix + tabName
    confRMSDsByClusMean = popResults['confRMSDsByClusMean']
    confRMSDsByClusSerr = popResults['confRMSDsByClusSerr']
    rows = len(confRMSDsByClusMean)
    cols = len(confRMSDsByClusMean[0])
    strVals = []
    for row in range(nConfs):
        vals = []
        for col in range(nMajorPlus1):
            mean = confRMSDsByClusMean[row][col]
            serr = confRMSDsByClusSerr[row][col]
            vals.append('{:4.2f} &plusmn; {:4.2f} &#8491'.format(mean, 1.96 * serr))
        strVals.append(vals)
    #print(strVals)
    # assemble dict of html components to go into the table formatting function
    html_ConfRMSDsByClus = dict()
    html_ConfRMSDsByClus['title'] = 'Cluster RMSD from each Starting Pose'
    html_ConfRMSDsByClus['subtitle'] = 'How close is each ligand cluster to each starting pose?'
    maxFieldChars = 12
    colorBy = 'column'
    style, body = HtmlTableType1(grdName, colNames, poseNames, strVals, colColors, colorBy, maxFieldChars)
    html_ConfRMSDsByClus['style'] = style
    html_ConfRMSDsByClus['body'] = body
    descr = 'This Table tells us how far each major cluster lies from each starting pose '
    descr += 'by computing the ensemble average RMS deviation '
    descr += '(in &#8491;) between them.'
    descr += 'The "Other" category includes minor clusters as well as Outliers from the clustering.'
    html_ConfRMSDsByClus['descr'] = descr


    # Modify the headers and text to center and have a decent max-width
    headermodsName = 'STMD_Clus_headermods'
    headermodsStyle = '\n.' + headermodsName + ' {\n  max-width: 800px;\n  margin: auto;\ntext-align: center;\n}\n\n'
    # print(headermodsStyle)
    htmlStyle = headermodsStyle
    # Omitting the html_ClusInfo Table because the data is already given elsewhere
    #htmlStyle += html_ClusInfo['style']
    htmlStyle += html_ClusFracByConf['style']
    htmlStyle += html_ConfRMSDsByClus['style']
    htmlBody = ''
    # Omitting the html_ClusInfo Table because the data is already given elsewhere
    #for html_dict in [html_ClusInfo, html_ClusFracByConf, html_ConfRMSDsByClus]:
    for html_dict in [html_ClusFracByConf, html_ConfRMSDsByClus]:
        htmlBody += HtmlWriteTableBody(html_dict, headermodsName)
    #print(htmlStyle)
    #print(htmlBody)

    return htmlStyle, htmlBody

