# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
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


# Acknowledgement
# The code has been adapted from the MIT licensed package
# https://github.com/choderalab/freeenergyframework
# developed by Hannah Bruce MacDonald from the Chodera Lab

import numpy as np
import sklearn.metrics
import scipy


def calc_RAE(y_true_sample, y_pred_sample):
    MAE = sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
    mean = np.mean(y_true_sample)
    MAD = np.sum([np.abs(mean - i) for i in y_true_sample]) / float(len(y_true_sample))
    return MAE / MAD


def calc_RRMSE(y_true_sample, y_pred_sample):
    mse = sklearn.metrics.mean_squared_error(y_true_sample, y_pred_sample)
    mean_exp = np.mean(y_true_sample)
    mds = np.sum([(mean_exp - i) ** 2 for i in y_true_sample]) / float(len(y_true_sample))
    rrmse = np.sqrt(mse / mds)

    return rrmse


def compute_statistic(y_true_sample, y_pred_sample, statistic):
    """Compute requested statistic.

    Parameters
    ----------
    y_true : ndarray with shape (N,)
        True values
    y_pred : ndarray with shape (N,)
        Predicted values
    statistic : str
        Statistic, one of ['RMSE', 'MUE', 'R2', 'rho']

    """

    if statistic == 'MAE':
        return sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
    elif statistic == 'RAE':
        return calc_RAE(y_true_sample, y_pred_sample)
    elif statistic == 'RMSE':
        return np.sqrt(sklearn.metrics.mean_squared_error(y_true_sample, y_pred_sample))
    elif statistic == 'RRMSE':
        return calc_RRMSE(y_true_sample, y_pred_sample)
    elif statistic == 'R2':
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_true_sample, y_pred_sample)
        return r_value**2
    elif statistic == 'RHO':
        return scipy.stats.pearsonr(y_true_sample, y_pred_sample)[0]
    elif statistic == 'KTAU':
        return scipy.stats.kendalltau(y_true_sample, y_pred_sample)[0]
    else:
        raise Exception("unknown statistic '{}'".format(statistic))


def bootstrap_statistic(y_true, y_pred, ci=0.95, statistic='RMSE', nbootstrap=1000, plot_type='dG'):

    """Compute mean and confidence intervals of specified statistic.

    Parameters
    ----------
    y_true : ndarray with shape (N,)
        True values
    y_pred : ndarray with shape (N,)
        Predicted values
    ci : float, optional, default=0.95
        Interval for CI
    statistic : str
        Statistic, one of ['RMSE', 'MUE', 'R2', 'rho']
    nbootstrap : int, optional, default=1000
        Number of bootstrap samples
    plot_type : str, optional, default='dG'
        'dG' or 'ddG'

    Returns
    -------
    rmse_stats : dict of floeat
        'mean' : mean RMSE
        'stderr' : standard error
        'low' : low end of CI
        'high' : high end of CI
    """
    def unique_differences(x):
        """Compute all unique differences"""
        N = len(x)
        return np.array([(x[i] - x[j]) for i in range(N) for j in range(N) if (i != j) ] )

    assert len(y_true) == len(y_pred)
    sample_size = len(y_true)
    s_n = np.zeros([nbootstrap], np.float64) # s_n[n] is the statistic computed for bootstrap sample n
    for replicate in range(nbootstrap):
        indices = np.random.choice(np.arange(sample_size), size=[sample_size])
        if plot_type == 'dG':
            y_true_sample, y_pred_sample = y_true[indices], y_pred[indices]
        elif plot_type == 'ddG':
            y_true_sample, y_pred_sample = unique_differences(y_true[indices]), unique_differences(y_pred[indices])
        s_n[replicate] = compute_statistic(y_true_sample, y_pred_sample, statistic)

    rmse_stats = dict()
    if plot_type == 'dG':
        rmse_stats['mle'] = compute_statistic(y_true, y_pred, statistic)
    elif plot_type == 'ddG':
        rmse_stats['mle'] = compute_statistic(unique_differences(y_true), unique_differences(y_pred), statistic)
    rmse_stats['stderr'] = np.std(s_n)
    rmse_stats['mean'] = np.mean(s_n)
    # TODO: Is there a canned method to do this?
    s_n = np.sort(s_n)
    low_frac = (1.0-ci)/2.0
    high_frac = 1.0 - low_frac
    rmse_stats['low'] = s_n[int(np.floor(nbootstrap*low_frac))]
    rmse_stats['high'] = s_n[int(np.ceil(nbootstrap*high_frac))]

    return rmse_stats

