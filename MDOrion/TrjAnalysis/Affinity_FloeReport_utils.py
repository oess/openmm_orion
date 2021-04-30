# (C) 2021 OpenEye Scientific Software Inc. All rights reserved.
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

import numpy as np
import math

import plotly
import plotly.graph_objs as go

from sklearn import linear_model
import MDOrion.TrjAnalysis.stats_utils as trjstats

_affnty_floe_report_header = '''<html>
<head>
<meta charset="utf-8" />
<style>

  .affinity-report-wrapper * {
    box-sizing: border-box;
    font-family: 'Helvetica Neue', Helvetica;
  }

  .affinity-report__row {
    width:100%;
    display:flex;
    align-items: center;
  }
  .affinity-report__half {
    width: 50%;
  }
  .affinity-report__60pct {
    width: 60%;
  }
  .affinity-report__40pct {
    width: 40%;
  }
  .affinity-report__column > * {
    width: 100%;
    height: auto;
  }

.AffinityModel_headermods {
  max-width: 800px;
  margin: auto;
text-align: center;
}

.AffinityModelTable-bold {
  font-weight: bold;
}

AffinityModelTable-item  {
  background-color: #fafafa;
  text-align: center;
  padding: 4px 0;
  font-size: 16px;
}

'''


_statsTableStyle = '''AffinityModelStatsTable-container {
  display: grid;
  grid-template-columns: auto auto auto auto;
  grid-template-rows: 1fr 1fr 1fr 1fr 1fr;
  grid-row-gap: 2px;
  grid-column-gap: 1px;
  padding: 2px;
  max-width: 810px;
  margin: auto;
  background-color: #a0a0a0;
}

'''


_RLModelTableStyle = '''LinearModelStatsTable-container {
  display: grid;
  grid-template-columns: auto auto auto;
  /* grid-template-rows: 1fr 1fr 1fr 1fr 1fr; */
  grid-row-gap: 2px;
  grid-column-gap: 1px;
  padding: 2px;
  max-width: 810px;
  margin: auto;
  background-color: #a0a0a0;
}'''


_affnty_style_footer = '''
</style>
</head>
<body>
'''

_AffinityTableSectionHeader = '  <div class="affinity-report__column affinity-report__half">\n\n'

_section_divider = '  <hr style="max-width: 1000px;">\n'

_html_trailer = '\n</body>\n</html>'


def GenerateGraphTablesHeaderHTML(metric_name):
    header = '''
    <div class="affinity-report-wrapper"> 
      <br><h2 style="text-align: center;">Predicted  {} vs. Experiment</h2>
      <div class="affinity-report__row">
      <div class="affinity-report__column affinity-report__half"
           style="height: 500px;">\n'''.format(metric_name)
    return header


def TrimPlotlyHTML( plotly_figure):
    # find start position for trim based on <body> string
    bodypos1 = plotly_figure.find('<body>', 0, 100)
    start_trim = bodypos1+6
    #print(len(plotly_figure),start_trim)
    #print(plotly_figure[:start_trim])
    #
    # find end position for trim based on </body> string
    bodypos2 = plotly_figure.find('</body>', len(plotly_figure)-100)
    end_trim = bodypos2-len(plotly_figure)
    #print(len(plotly_figure),bodypos2, end_trim)
    #print(plotly_figure[end_trim:])
    #
    # return substring in between these positions
    return plotly_figure[start_trim:end_trim]


def PlotlyLinearModelOfX(xdata, ydata, slope, intercept,lineColor='red',label=None):
    # Interpolate range of x
    pltx = np.linspace(math.floor(np.min(xdata)), math.ceil(np.max(xdata)), 100, endpoint=True)
    plty = trjstats.xFromYLinearModel(pltx,slope,intercept)

    #
    # exclude interpolated points that fall outside ydata range
    TSLine_x = []
    TSLine_y = []
    for x,y in zip(pltx,plty):
        if y>=math.floor(np.min(ydata)) and y<=math.ceil(np.max(ydata)):
            TSLine_x.append(x)
            TSLine_y.append(y)
    #print(len(TSLine_x),len(TSLine_y))
    #
    if label is not None:
        hovtext = '{}<br>x = {:.3f} * y + {:.3f}'.format(label,slope,intercept)
        gname = label
        showlegd=True
    else:
        hovtext = 'x = {:.3f} * y + {:.3f}'.format(slope,intercept)
        gname = ''
        showlegd=False
    # generate the dynamically labeled interpolation line
    interpolation = go.Scatter(
        name=gname,
        showlegend=showlegd,
        x=TSLine_x,
        y=TSLine_y,
        hovertext=hovtext,
        hoverinfo="text",
        mode="lines",
        line=dict(color=lineColor, width=2))
    #
    return interpolation


def MakeAffinityModelFigure(figure, xdata,err_x,ydata,err_y,hover_text,hov_label):
    # Generate the basic graph of the data with hover text
    trace0 = go.Scatter(
        x = xdata, y = ydata,
        text=hover_text,
        hovertemplate=hov_label,
        name='',
        mode='markers',
        marker=dict(color='purple', size=8),
        xaxis='x1',
        yaxis='y1',
        showlegend=False,
        error_x=dict(
            type='data',
            array=err_x,
            color='purple',
            visible=True),
        error_y=dict(
            type='data',
            array=err_y,
            color='purple',
            visible=True),
    )
    figure.add_trace(trace0)
    #
    return figure


def AffinityModelFigureAddLayout(figure,xaxlabel,yaxlabel,square_plot_axes=None):
# Generate axes with labels
    xaxis1 = dict(
            zeroline=True,
            showgrid=True,
            tickfont=dict(
                size=18,
            ),
            constrain='domain',
            title=xaxlabel,
            titlefont=dict(
                size=20
            )
        )
    # square_plot_axes defines a range for the x and y axes when a square plot is desired
    if square_plot_axes is not None:
        xaxis1['range'] = square_plot_axes

    yaxis1 = dict(
            zeroline=True,
            showgrid=True,
            tickfont=dict(
                size=18,
            ),
            title=yaxlabel,
            titlefont=dict(
                size=20
            )
        )
    if square_plot_axes is not None:
        yaxis1['range'] = square_plot_axes
        yaxis1['scaleanchor'] = "x1"
        #yaxis1['scaleratio'] = 1
    #
    # Refine the graph layout
    figure.update_layout(
        xaxis1=xaxis1,
        yaxis1=yaxis1,
        margin=dict(
            l=120,
            r=30,
            b=50,
            t=30,
            pad=4),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01)
        )
    #
    return figure


def FreeEnergyBasePlot(xrange, kcal_err=1):

    x_plt = np.linspace(xrange[0], xrange[1], 100, endpoint=True)

    x_plus_1kcal = x_plt + kcal_err
    x_minus_1kcal = x_plt - kcal_err

    alpha_color = 'rgba(180, 180, 180, 0.4)'

    figure = go.Figure()

    # Theoretical line minus 1kcal
    figure.add_trace(go.Scatter(
        x=x_plt,
        y=x_minus_1kcal,
        fill=None,
        hovertext="",
        hoverinfo="text",
        mode="lines",
        showlegend=False,
        line=dict(color=alpha_color, width=1))
    )

    # Theoretical line plus 1kcal
    figure.add_trace(go.Scatter(
        x=x_plt,
        y=x_plus_1kcal,
        fill='tonexty',
        fillcolor=alpha_color,
        hovertext="",
        hoverinfo="text",
        mode="lines",
        showlegend=False,
        line=dict(color=alpha_color, width=1))
    )

    # Theoretical line
    figure.add_trace(go.Scatter(
        x=x_plt,
        y=x_plt,
        hovertext="1:1",
        hoverinfo="text",
        mode="lines",
        showlegend=False,
        line=dict(color='black', width=2, dash='dash'))
    )

    return figure


def HorizLinesDGPredicted(figure, y, err_y, hover_text, xrange):
    x_plt = np.linspace(xrange[0], xrange[1], 20, endpoint=True)
    nx = len(x_plt)

    y_plus_err = x_plt + 1
    y_minus_err = x_plt - 1

    hov_label = "<b>Ligand %{text}</b><br>" + "Pred \u0394G: %{y:.2f}<br>"
    name_plot = '\u0394G'

    alpha_color = 'rgba(243, 216, 216, 0.5)'
    # alpha_color = 'BlanchedAlmond'

    for yi, yi_err, ligname in zip(y, err_y, hover_text):
        y_plt = [yi] * nx
        y_plus = [yi + yi_err] * nx
        y_minus = [yi - yi_err] * nx
        htext = [ligname] * nx

        # Theoretical line minus 1kcal
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=y_minus,
            fill=None,
            hovertext="",
            hoverinfo="text",
            mode="lines",
            showlegend=False,
            line=dict(color=alpha_color, width=1))
        )

        # Theoretical line plus 1kcal
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=y_plus,
            fill='tonexty',
            fillcolor=alpha_color,
            hovertext="",
            hoverinfo="text",
            mode="lines",
            showlegend=False,
            line=dict(color=alpha_color, width=1))
        )

        # Theoretical line
        figure.add_trace(go.Scatter(
            x=x_plt,
            y=y_plt,
            text=htext,
            hovertemplate=hov_label,
            name=name_plot,
            mode="lines",
            showlegend=False,
            line=dict(color='IndianRed', width=2))
        )

    return figure


def FindAxisRangeForSquarePlot(x, xerr, y, yerr):
    # Axis range:
    conc_axes = np.concatenate((x, y))
    conc_daxes = np.concatenate((xerr, yerr))

    # Correct for the error bars
    min_range_axes = np.min(conc_axes - conc_daxes)
    max_range_axes = np.max(conc_axes + conc_daxes)
    # print('min, max range_axes', min_range_axes, max_range_axes)

    padding = (max_range_axes - min_range_axes) * 0.05
    # print('padding', padding)

    min_range_axes -= padding
    max_range_axes += padding
    # print('padded min, max range_axes', min_range_axes, max_range_axes)

    return [min_range_axes, max_range_axes]


def CombinePredExptDDGs(lig_expt,edge_pred):
    lig_names = set(lig_expt.keys())
    combined = dict()
    for edge_name in edge_pred.keys():
        lig_name_state_A = edge_name.split()[0]
        lig_name_state_B = edge_name.split()[2]
        combined[edge_name] = [edge_pred[edge_name][0], edge_pred[edge_name][1]]
        if lig_name_state_A in lig_names and lig_name_state_B in lig_names:
            # Experimental relative binding affinity in kJ/mol
            DDG_A_to_B_exp = lig_expt[lig_name_state_B][0] - lig_expt[lig_name_state_A][0]
            # Experimental relative binding affinity error
            ddG_A_to_B_exp = math.sqrt(lig_expt[lig_name_state_B][1]**2 + lig_expt[lig_name_state_A][1]**2)
            combined[edge_name] = combined[edge_name] + [DDG_A_to_B_exp, ddG_A_to_B_exp]
        else:
            combined[edge_name] = combined[edge_name] + [None, None]
    #
    return combined


def CombineAndOffsetPredExptDGs(lig_expt,affinity_pred):
    lig_names = set(lig_expt.keys())
    combined = dict()
    for lig_name in affinity_pred.keys():
        combined[lig_name] = affinity_pred[lig_name]
        if lig_name in lig_names:
            # Append experimental binding affinity & error
            DG_A_to_B_exp = lig_expt[lig_name][0]
            # Experimental binding affinity error
            dG_A_to_B_exp = lig_expt[lig_name][1]
            combined[lig_name] = combined[lig_name] + [lig_expt[lig_name][0], lig_expt[lig_name][1]]
        else:
            combined[lig_name] = combined[lig_name] + [None, None]
    #
    # select data with both predicted and experimental data to figure out DG offset
    exp_DG = []
    pred_DG_with_expt = []
    for lig_name in combined.keys():
        lig = combined[lig_name]
        if lig[2] is not None:
            pred_DG_with_expt.append(lig[0])
            exp_DG.append(lig[2])
    #print(exp_DG)
    #print(pred_DG_with_expt)
    #
    # Offset Predicted DGs to have the same mean as that of the experimental DGs
    exp_DG_mean = sum(exp_DG)/len(exp_DG)
    pred_DG_arr = np.array(pred_DG_with_expt)
    offset = exp_DG_mean - pred_DG_arr.mean()
    #print(exp_DG_mean, offset)
    #
    for lig_name in combined.keys():
        combined[lig_name][0] += offset
    #
    return combined


def GeneratePredExptBothDatavecs(pred_expt_dic, units='kcal/mol', symmetrize=False):
    # Extract separate datavecs for items having both predicted and experimental data

    if units == 'kcal/mol':
        conv_factor = 1 / 4.184
    else:
        conv_factor = 1.0

    # Extract data only for items which have both predicted and expt data
    pred = []
    dpred = []
    expt = []
    dexpt = []
    #
    for name, item in pred_expt_dic.items():
        if item[2] is not None:
            pred.append(conv_factor * item[0])
            dpred.append(conv_factor * item[1])
            expt.append(conv_factor * item[2])
            dexpt.append(conv_factor * item[3])
            if symmetrize:
                pred.append(conv_factor * -item[0])
                dpred.append(conv_factor * item[1])
                expt.append(conv_factor * -item[2])
                dexpt.append(conv_factor * item[3])
    #
    return pred, dpred, expt, dexpt


def Generate_DDG_Plot(edge_pred_expt_dic,units='kcal/mol',symmetrize=True):
    #
    if units == 'kcal/mol':
        conv_factor = 1/4.184
    else:
        conv_factor = 1.0
    #
    #Aggregate plot data only for edges which have both predicted and expt data
    pred_DDG = []
    dpred_DDG = []
    exp_DDG = []
    dexp_DDG = []
    hov_edge = []
    #
    for edge_name in edge_pred_expt_dic.keys():
        lig_name_state_A = edge_name.split()[0]
        lig_name_state_B = edge_name.split()[2]
        edge = edge_pred_expt_dic[edge_name]
        if edge[2] is not None:
            pred_DDG.append(conv_factor * edge[0])
            dpred_DDG.append(conv_factor * edge[1])
            exp_DDG.append(conv_factor * edge[2])
            dexp_DDG.append(conv_factor * edge[3])
            hov_edge.append('{} to {}'.format(lig_name_state_A, lig_name_state_B))
            if symmetrize:
                pred_DDG.append(conv_factor * -edge[0])
                dpred_DDG.append(conv_factor * edge[1])
                exp_DDG.append(conv_factor * -edge[2])
                dexp_DDG.append(conv_factor * edge[3])
                hov_edge.append('{} to {}'.format(lig_name_state_B, lig_name_state_A))
    #
    range_axis_DDG = FindAxisRangeForSquarePlot(exp_DDG, dexp_DDG, pred_DDG, dpred_DDG)
    #print(range_axis_DDG)
    #
    fig_DDG_base = FreeEnergyBasePlot(range_axis_DDG)
    #
    hov_label = "<b>%{text}</b><br>" + "Expt \u0394\u0394G: %{x:.2f}<br>" + "Pred \u0394\u0394G: %{y:.2f}<br>"
    fig_DDG_plt = MakeAffinityModelFigure(fig_DDG_base, exp_DDG, dexp_DDG, pred_DDG, dpred_DDG,
                                          hov_edge, hov_label)
    #
    xlabel = 'Experimental \u0394\u0394G ({})'.format(units)
    ylabel = 'Predicted \u0394\u0394G ({})'.format(units)
    fig_DDG = AffinityModelFigureAddLayout(fig_DDG_plt,xlabel,ylabel,range_axis_DDG)
    #
    return fig_DDG


def Generate_DG_Plot(lig_pred_expt_dic,units='kcal/mol'):
    # Set constants
    display_units = units
    if units == 'kcal/mol':
        conv_factor = 1/4.184
    else:
        conv_factor = 1.0
    # Aggregate plot data only for edges which have both predicted and expt data
    pred_DG = []
    dpred_DG = []
    exp_DG = []
    dexp_DG = []
    hov_ligs = []
    #
    for lig_name in lig_pred_expt_dic.keys():
        lig = lig_pred_expt_dic[lig_name]
        #print('lig', lig)
        if lig[2] is not None:
            pred_DG.append(conv_factor * lig[0])
            dpred_DG.append(conv_factor * lig[1])
            exp_DG.append(conv_factor * lig[2])
            dexp_DG.append(conv_factor * lig[3])
            hov_ligs.append(lig_name)
    #
    range_axis_DG = FindAxisRangeForSquarePlot(exp_DG, dexp_DG, pred_DG, dpred_DG)
    #print(range_axis_DG)
    #
    fig_DG_base = FreeEnergyBasePlot(range_axis_DG)
    #
    hov_label = "<b>%{text}</b><br>" + "Expt \u0394G: %{x:.2f}<br>" + "Pred \u0394G: %{y:.2f}<br>"
    fig_DG_pred = MakeAffinityModelFigure(fig_DG_base, exp_DG, dexp_DG, pred_DG, dpred_DG,
                                          hov_ligs, hov_label)
    #
    xlabel = 'Experimental \u0394G {}'.format(display_units)
    ylabel = 'Predicted \u0394G {}'.format(display_units)
    fig_DG = AffinityModelFigureAddLayout(fig_DG_pred, xlabel, ylabel, range_axis_DG)
    #
    return fig_DG


def GeneratePredDGTable(lig_pred_expt, units='kcal/mol', glabel='\u0394G'):
    if units == 'kcal/mol':
        conv_factor = 1/4.184
    else:
        conv_factor = 1.0
        #
    tabName = 'AffinityModelTable'
    if glabel=='\u0394\u0394G':
        headers = ['RBFE Edge', 'Predicted<sup>a</sup>', 'Experiment<sup>a</sup>']
    else:
        headers = ['Ligand Name', 'Predicted<sup>a,b</sup>', 'Experiment<sup>b</sup>']
    #
    # table header
    table = '  <br><h2 class="AffinityModel_headermods">Binding {}: Predicted and Experimental</h2>'.format(glabel)
    # table body
    strtBold = '  <'+tabName+'-item class="AffinityModelTable-bold">'
    strt = '  <'+tabName+'-item>'
    end = '</'+tabName+'-item>\n'
    table += '\n<LinearModelStatsTable-container>\n'
    for header in headers:
        table += strtBold+header+end
    #
    for lig in sorted(lig_pred_expt.items(),key=lambda x:x[1]):
        table += strtBold+'{}'.format(lig[0])+end
        dg = conv_factor*lig[1][0]
        ddg = conv_factor*lig[1][1]
        table += strt+'{:.1f} &plusmn; {:.1f}'.format(dg, ddg)+end
        if lig[1][2] is None:
            table += strt+'N/A <sup>c</sup>'+end
        else:
            dg = conv_factor*lig[1][2]
            ddg = conv_factor*lig[1][3]
            table += strt+'{:.1f} &plusmn; {:.1f}'.format(dg, ddg)+end
    #
    table+= '</LinearModelStatsTable-container>\n'
    #
    # table footer
    if glabel=='\u0394\u0394G':
        table+= '''  <p class=\"AffinityModel_headermods\"> 
        <sup>a</sup>Values in {}.
        <sup>c</sup>No experimental data provided.<br>\n\n'''.format(units)
    else:
        table+= '''  <p class=\"AffinityModel_headermods\">
        <sup>a</sup>Maximum Likelihood Estimation using MBAR; 
        Offset using the mean of the experimental data as a reference. 
        <sup>b</sup>Values in {}.
        <sup>c</sup>No experimental data provided.<br>\n\n'''.format(units)
    #
    return table


def GenerateCorrTable(statsResult, xlabel, ylabel, boot_nreps, units='kcal/mol', showMAE=False):
    tabName = 'AffinityModelTable'
    colNames = ['Metric', 'Value [5%,95%]<sup>a</sup>', 'StdErr<sup>a</sup>', '<i>p</i>-value']
    metrics = [ 'Pearsons r-sq', 'Kendalls tau', 'MAE', 'RMAE']
    metricHTML = {'Pearsons r-sq':"Pearson's <i>r</i>&sup2;",
                  'Kendalls tau':"Kendall's &tau;",
                  'MAE': 'MAE<sup>b</sup>',
                  'RMAE': 'RMAE<sup>c</sup>'
                 }
    #
    # table header
    table = '  <br><h3 class="AffinityModel_headermods">Correlations Between {} and {}</h3>'.format(xlabel,ylabel)
    # table body
    strtBold = '  <'+tabName+'-item class="AffinityModelTable-bold">'
    strt = '  <'+tabName+'-item>'
    end = '</'+tabName+'-item>\n'
    table += '\n<AffinityModelStatsTable-container>\n'
    for field in colNames:
        table += strtBold+field+end
    #
    #for name, rowVals in zip(metrics, tabVals):
    for metric in metrics:
        if metric not in metricHTML.keys():
            continue
        if not showMAE and (metric=='MAE' or metric=='RMAE'):
            continue
        table += strtBold+metricHTML[metric]+end
        table += (strt+'{:.3f}'.format(statsResult[metric]['value'])+' ['
                    +'{:.3f}'.format(statsResult[metric]['bootconf05'])+', '
                    +'{:.3f}'.format(statsResult[metric]['bootconf95'])+']'+end)
        table += strt+'{:.3f}'.format(statsResult[metric]['bootstd'])+end
        pval = statsResult[metric]['p-value']
        if pval:
            table += strt+'{:.3e}'.format(statsResult[metric]['p-value'])+end
        else:
            table += strt+'n/a'+end
    #
    table += '</AffinityModelStatsTable-container>\n'
    #
    # table footer
    table += '''  <p class=\"AffinityModel_headermods\">
    <sup>a</sup> 5% to 95% Confidence Intervals and Standard Error estimated using bootstrapping ({} samples).
    <sup>b</sup>Mean Absolute Error in {}. 
    <sup>c</sup>MAE divided by the Mean Absolute Deviation of {}.</p><br>\n\n'''.format(boot_nreps, units, xlabel)
    #
    return table

def GenerateRLModelTable(RLModelInfo, xlabel, ylabel, units='kcal/mol'):
    tabName = 'AffinityModelTable'
    models = ['Theil-Sen', 'Huber']
    metricHTML = {'slope':'Slope',
                  'intcp':'Intercept',
                  'MAE':"MAE<sup>a</sup>",
                  'RMAE':"Relative MAE<sup>b</sup>"}
    #
    # table header
    table = '  <br><h3 class="AffinityModel_headermods">Robust Linear Models of {}</h3>'.format(xlabel)
    # table body
    strtBold = '  <'+tabName+'-item class="AffinityModelTable-bold">'
    strt = '  <'+tabName+'-item>'
    end = '</'+tabName+'-item>\n'
    table += '\n<LinearModelStatsTable-container>\n'
    table += strtBold+'Metric'+end
    for model in models:
        table += strtBold+model+end
    #for name, rowVals in zip(metrics, tabVals):
    #for metric in metrics:
    for metric in metricHTML:
        table += strtBold+metricHTML[metric]+end
        table += strt+'{:.3f}'.format(RLModelInfo['Theil-Sen '+metric])+end
        table += strt+'{:.3f}'.format(RLModelInfo['Huber '+metric])+end
    table += '</LinearModelStatsTable-container>\n'
    #
    # table footer
    table += '''  <p class=\"AffinityModel_headermods\">
    <sup>a</sup>Mean Absolute Error in {}. 
    <sup>b</sup>MAE divided by the Absolute Average Deviation of {}.</p><br>\n\n'''.format(units, xlabel)
    #
    return table


def GenerateAffinityGraphTablesHTML(pred_expt_dic, glabel,
                                    units='kcal/mol',
                                    symmetrize=False,
                                    nreps=1000,
                                    showMAE=True):

    # Extract data only for ligands which have both predicted and expt data
    pred, dpred, expt, dexpt = GeneratePredExptBothDatavecs(pred_expt_dic, units=units, symmetrize=symmetrize)

    # Calculate Correlation Stats, raw MAE, and raw Relative MAE
    statsResult = trjstats.PearsonKendallMAEStats(expt, pred)

    # Bootstrap confidence intervals and std error
    bootResult = trjstats.BootStats(expt, pred, nreps=nreps)
    for metric in statsResult.keys():
        for bootmetric in bootResult[metric].keys():
            statsResult[metric][bootmetric] = bootResult[metric][bootmetric]

    # Generate html table for correlation/MAE table
    xlabel = 'Experimental {}'.format(glabel)
    ylabel = 'Predicted {}'.format(glabel)
    html_statsTable = GenerateCorrTable(statsResult, xlabel, ylabel, nreps, units=units, showMAE=showMAE)

    # Robust linear estimators of deltaG
    lModResult = trjstats.GenerateHuberTheilSenRLModels(pred, expt)

    # Generate html table for Robust Linear (RL) models
    html_RLModelTable = GenerateRLModelTable(lModResult, xlabel, ylabel, units=units)

    # Generate the DG plot and add lines for RL models
    if glabel=='\u0394\u0394G':
        fig = Generate_DDG_Plot(pred_expt_dic,units=units,symmetrize=symmetrize)
    else:
        fig = Generate_DG_Plot(pred_expt_dic, units=units)

    range_axis_DG = FindAxisRangeForSquarePlot(pred, dpred, expt, dexpt)

    trace_TheilSen = PlotlyLinearModelOfX(range_axis_DG, range_axis_DG,
                                                lModResult['Theil-Sen slope'],
                                                lModResult['Theil-Sen intcp'],
                                                label='Theil-Sen Estimator', lineColor='red')

    trace_Huber = PlotlyLinearModelOfX(range_axis_DG, range_axis_DG,
                                             lModResult['Huber slope'],
                                             lModResult['Huber intcp'],
                                             label='Huber Estimator', lineColor='orange')

    fig.add_trace(trace_TheilSen)
    fig.add_trace(trace_Huber)
    fig_fixLegend = fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.20))

    fig_html_base = fig_fixLegend.to_html()
    html_fig = TrimPlotlyHTML(fig_html_base)

    return html_statsTable, html_RLModelTable, html_fig