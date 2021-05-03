# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2017 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------


import sys

from scipy.integrate import simps

from scipy.optimize import fmin

from scipy import stats

import numpy as np


def _longest_dgdl_file(lst):
    """Takes a list of dgdl.xvg files and returns the id (starting from 0) of
    of the longest file.
    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.
    Returns
    -------
    ind : int
        index of the longest file
    """

    ind = 0
    maxlen = -9999
    for idx, f in enumerate(lst[0:]):
        fp = open(f)
        lines = fp.readlines()
        if len(lines) > maxlen:
            maxlen = len(lines)
            ind = idx

    return ind


def parse_dgdl_files(lst, lambda0=0, invert_values=False):
    """Takes a list of dgdl.xvg files and returns the integrated work values
    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.
    Returns
    -------
    w : array
        array of work values.
    lst : list
        sorted list of input dgdl.avg files corresponding to the work values
        in w.
    """

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    # identify file with the most entries
    imax = _longest_dgdl_file( lst )

    _check_dgdl(lst[imax], lambda0)
    first_w, ndata = integrate_dgdl(lst[imax], lambda0=lambda0,
                                    invert_values=invert_values)
    w_list = [first_w]
    for idx, f in enumerate(lst[0:]):
        if idx == imax:  # already processed
            continue
        sys.stdout.write('\r    Reading %s' % f)
        sys.stdout.flush()

        w, _ = integrate_dgdl(f, ndata=ndata, lambda0=lambda0,
                              invert_values=invert_values)
        if w is not None:
            w_list.append(w)

    print('\n')

    return w_list


def integrate_dgdl(fn, ndata=-1, lambda0=0, invert_values=False):
    """Integrates the data in a dgdl.xvg file.
    Parameters
    ----------
    fn : str
        the inpur dgdl.xvg file from Gromacs.
    ndata : int, optional
        number of datapoints in file. If -1, then ??? default is -1.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.
    Returns
    -------
    integr : float
        result of the integration performed using Simpson's rule.
    ndata : int
        number of data points in the input file.
    """

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    lines = open(fn).readlines()
    if not lines:
        return None, None

    # extract dgdl datapoints into r
    # TODO: we removed the check for file integrity. We could have an
    # optional files integrity check before calling this integration func

    lines = [l for l in lines if l[0] not in '#@&']
    try:
        r = [float(x.split()[1]) for x in lines]
    except:
        r = "incomplete_file"
        ndata = 1

    if ndata != -1 and len(r) != ndata:
        try:
            if "incomplete_file" in r:
                print(' !! Skipping %s (incomplete file, probably simulation crashed)\n' % fn)
            else:
                print(' !! Skipping %s ( read %d data points, should be %d )' % (fn, len(r), ndata))
        except:
            print(' !! Skipping %s ' % (fn))
        return None, None
    # convert time to lambda
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1

    # arrays for the integration
    # --------------------------
    # array of lambda values
    x = [lambda0+i*dlambda for i, dgdl in enumerate(r)]
    # array of dgdl
    y = r

    if lambda0 == 1:
        x.reverse()
        y.reverse()

    if invert_values is True:
        integr = simps(y, x) * (-1)
        return integr, ndata
    else:
        integr = simps(y, x)
        return integr, ndata


def _check_dgdl(fn, lambda0):
    """Prints some info about a dgdl.xvg file."""
    l = open(fn).readlines()
    if not l:
        return None
    r = []
    for line in l:
        if line[0] not in '#@&':
            r.append([float(x) for x in line.split()])
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1

    print('    # data points: %d' % ndata)
    print('    Length of trajectory: %8.3f ps' % r[-1][0])
    print('    Delta lambda: %8.5f' % dlambda)


# Constants
kb = 0.00831447215   # kJ/(K*mol)


class BAR(object):
    """Bennett acceptance ratio (BAR).
    Description...
    Parameters
    ----------
    Examples
    --------
    """

    def __init__(self, wf, wr, T, nboots=0, nblocks=1):
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots
        self.nblocks = nblocks

        self.nf = len(wf)
        self.nr = len(wr)
        self.beta = 1./(kb*self.T)
        self.M = kb * self.T * np.log(float(self.nf) / float(self.nr))

        # Calculate all BAR properties available
        self.dg = self.calc_dg(self.wf, self.wr, self.T)
        self.err = self.calc_err(self.dg, self.wf, self.wr, self.T)
        if nboots > 0:
            self.err_boot = self.calc_err_boot(self.wf, self.wr, nboots,
                                               self.T)
        self.conv = self.calc_conv(self.dg, self.wf, self.wr, self.T)
        if nboots > 0:
            self.conv_err_boot = self.calc_conv_err_boot(self.dg, self.wf,
                                                         self.wr, nboots,
                                                         self.T)
        if nblocks > 1:
            self.err_blocks = self.calc_err_blocks(self.wf, self.wr, nblocks,
                                                   self.T)

    @staticmethod
    def calc_dg(wf, wr, T):
        """Estimates and returns the free energy difference.
        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        Returns
        ----------
        dg : float
            the BAR free energy estimate.
        """

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * np.log(nf/nr)

        def func(x, wf, wr):
            sf = 0
            for v in wf:
                sf += 1./(1+np.exp(beta*(M+v-x)))

            sr = 0
            for v in wr:
                sr += 1./(1+np.exp(-beta*(M+v-x)))

            r = sf-sr
            return r**2

        avA = np.average(wf)
        avB = np.average(wr)
        x0 = (avA+avB)/2.
        dg = fmin(func, x0=x0, args=(wf, wr), disp=0)

        return float(dg)

    @staticmethod
    def calc_err(dg, wf, wr, T):
        """Calculates the analytical error estimate.
        Parameters
        ----------
        dg : float
            the BAR free energy estimate
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        """

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * np.log(nf/nr)

        err = 0
        for v in wf:
            err += 1./(2+2*np.cosh(beta * (M+v-dg)))
        for v in wr:
            err += 1./(2+2*np.cosh(beta * (M+v-dg)))
        N = nf + nr
        err /= float(N)
        tot = 1/(beta**2*N)*(1./err-(N/nf + N/nr))

        err = float(np.sqrt(tot))
        return err

    @staticmethod
    def calc_err_boot(wf, wr, nboots, T):
        """Calculates the error by bootstrapping.
        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        nboots: int
            number of bootstrap samples.
        """

        nf = len(wf)
        nr = len(wr)
        dg_boots = []
        for k in range(nboots):
            # sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
            #                  % (k+1, nboots))
            # sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            dg_boot = BAR.calc_dg(bootA, bootB, T)
            dg_boots.append(dg_boot)

        print("Bootstrap iterations Stderr {} completed\n".format(nboots))

        err_boot = np.std(dg_boots)

        return err_boot

    @staticmethod
    def calc_err_blocks(wf, wr, nblocks, T):
        """Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.
        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        """

        dg_blocks = []
        # loosely split the arrays
        wf_split = np.array_split(wf, nblocks)
        wr_split = np.array_split(wr, nblocks)

        # calculate all dg
        for wf_block, wr_block in zip(wf_split, wr_split):
            dg_block = BAR.calc_dg(wf_block, wr_block, T)
            dg_blocks.append(dg_block)

        # get std err
        err_blocks = stats.sem(dg_blocks, ddof=1)

        return err_blocks

    @staticmethod
    def calc_conv(dg, wf, wr, T):
        """Evaluates BAR convergence as described in Hahn & Then, Phys Rev E
        (2010), 81, 041117. Returns a value between -1 and 1: the closer this
        value to zero the better the BAR convergence.
        Parameters
        ----------
        dg : float
            the BAR free energy estimate
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        """

        wf = np.array(wf)
        wr = np.array(wr)

        beta = 1./(kb*T)
        nf = len(wf)
        nr = len(wr)
        N = float(nf + nr)

        ratio_alpha = float(nf)/N
        ratio_beta = float(nr)/N
        bf = 1.0/(ratio_beta + ratio_alpha * np.exp(beta*(wf-dg)))
        tf = 1.0/(ratio_alpha + ratio_beta * np.exp(beta*(-wr+dg)))
        Ua = (np.mean(tf) + np.mean(bf))/2.0
        Ua2 = (ratio_alpha * np.mean(np.power(tf, 2)) +
               ratio_beta * np.mean(np.power(bf, 2)))
        conv = (Ua-Ua2)/Ua
        return conv

    @staticmethod
    def calc_conv_err_boot(dg, wf, wr, nboots, T):
        nf = len(wf)
        nr = len(wr)
        conv_boots = []
        for k in range(nboots):
            # sys.stdout.write('\r  Bootstrap (Conv): '
            #                  'iteration %s/%s' % (k+1, nboots))
            # sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            conv_boot = BAR.calc_conv(dg, bootA, bootB, T)
            conv_boots.append(conv_boot)

        print("Bootstrap iterations Conv {} completed\n".format(nboots))

        err = np.std(conv_boots)
        return err
