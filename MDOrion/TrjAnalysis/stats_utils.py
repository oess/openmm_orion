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

from scipy import stats
from sklearn import linear_model

def PearsonKendallMAEStats(xdata,ydata):
    corrResult = {}
    corrResult['Pearsons r-sq'] = {}
    corrResult['Kendalls tau'] = {}
    corrResult['MAE'] = {}
    corrResult['RMAE'] = {}
    prsnr0,prsnp0 = stats.pearsonr(xdata,ydata)
    corrResult['Pearsons r-sq']['value'] = prsnr0*prsnr0
    corrResult['Pearsons r-sq']['p-value'] = prsnp0
    #print('Pearsons r, r-sq, and p:', prsnr0,prsnr0*prsnr0,prsnp0)

    kndlltau0,kndllp0 = stats.kendalltau(xdata,ydata)
    corrResult['Kendalls tau']['value'] = kndlltau0
    corrResult['Kendalls tau']['p-value'] = kndllp0
    #print('Kendalls tau and p:', kndlltau0,kndllp0)

    # bootstrap RBFE MAE and RMAE
    MAE,RMAE = MeanAbsErrAndRelMAE(xdata,ydata)
    corrResult['MAE']['value'] = MAE
    corrResult['MAE']['p-value'] = None
    corrResult['RMAE']['value'] = RMAE
    corrResult['RMAE']['p-value'] = None

    return corrResult


def BootStats(xdata,ydata,nreps=2000):
    metrics = ['Pearsons r-sq', 'Kendalls tau', 'MAE', 'RMAE']
    boots = {}
    for metric in metrics:
        boots[metric] = []

    dataIdxs = list(range(len(xdata)))
    for rep in range(nreps):
        # bootstrap the indices of elements to pull from xdata and ydata
        bootIdx = np.random.choice(dataIdxs,size=len(dataIdxs))
        # skip this rep if all indices are identical: stats undefined
        if len(set(bootIdx))<2:
            continue
        # pull the bootstrapped indices from xdata and ydata
        xboot = [ xdata[i] for i in bootIdx]
        yboot = [ ydata[i] for i in bootIdx]
        # Claculate stats for the bootstrapped sample
        r, p = stats.pearsonr(xboot,yboot)
        boots['Pearsons r-sq'].append(r*r)
        tau,ptau = stats.kendalltau(xboot,yboot)
        boots['Kendalls tau'].append(tau)
        # bootstrap RBFE MAE and RMAE
        MAE,RMAE = MeanAbsErrAndRelMAE(xboot,yboot)
        boots['MAE'].append(MAE)
        boots['RMAE'].append(RMAE)

    bootResult = {}
    for metric in metrics:
        bootResult[metric] = {}
    for metric in metrics:
        #print()
        statsDict = StatsBootStats(boots[metric])
        #print(metric,statsDict)
        for bootmetric in statsDict.keys():
            bootResult[metric][bootmetric] = statsDict[bootmetric]

    return bootResult


def StatsBootStats(bootlist):
    arr = np.array(bootlist)
    stats = {}
    stats['bootmean'] = arr.mean()
    stats['bootmed'] = np.median(arr)
    stats['bootstd'] = arr.std()
    arr.sort()
    #print(len(arr),int(.05*len(arr)),int(.95*len(arr)))
    stats['bootconf05'] = arr[int(.05*len(arr))]
    stats['bootconf95'] = arr[int(.95*len(arr))]
    #
    return stats


def MeanAbsErrAndRelMAE(data,model):
    if len(data)!=len(model):
        return None
    vdata = np.array(data)
    vmodel = np.array(model)
    mad = np.abs(vdata-vdata.mean()).sum()/len(vdata)
    mae = np.abs(vmodel-vdata).sum()/len(vmodel)
    return mae, mae/mad


def xFromYLinearModel(y, slope, intercept):
    yarr = np.array(y)
    xarr = (yarr-intercept)/slope
    return xarr.tolist()


def yFromXLinearModel(x, slope, intercept):
    xarr = np.array(x)
    ypred = slope*xarr+intercept
    return ypred.tolist()


def RobustLinearModelWithStats(var_indep, var_dep, results, regressor, regressor_label=''):
    # add a space to the regressor_label
    if regressor_label != '':
        regressor_label += ' '
    # reformat independent variable as 1 dimensional feature vectors
    features = [[i] for i in var_indep]
    #print(features)
    #
    # fit robust linear model and add coefficients to results dictionary
    model = regressor.fit(features,var_dep)
    slope = model.coef_[0]
    intercept = model.intercept_
    results[regressor_label+'slope'] = slope
    results[regressor_label+'intcp'] = intercept
    #
    # compute MAE and RMAE and add results to results dictionary
    pred = yFromXLinearModel(var_indep, slope, intercept)
    MAE,RMAE = MeanAbsErrAndRelMAE(var_dep, pred)
    results[regressor_label+'MAE'] = MAE
    results[regressor_label+'RMAE'] = RMAE
    #
    return True


def GenerateHuberTheilSenRLModels(indep_var, dep_var):
    RLModResult = {}
    # Huber regressor
    hregr_label = 'Huber'
    hregr = linear_model.HuberRegressor()
    RobustLinearModelWithStats(indep_var, dep_var, RLModResult, hregr, hregr_label)
    #
    # Theil-Sen regressor
    tregr_label = 'Theil-Sen'
    tregr = linear_model.TheilSenRegressor()
    RobustLinearModelWithStats(indep_var, dep_var, RLModResult, tregr, tregr_label)
    #
    return RLModResult


