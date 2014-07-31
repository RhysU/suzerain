#!/usr/bin/env python
"""Usage: perfect.py [OPTIONS...] H5FILE...
Produces TBD

Options:
    -h            Display this help message and exit
    -o OUTSUFFIX  Save the output file PLOT.OUTSUFFIX instead of displaying

File H5FILE should have been produced by the perfect_summary application.
"""
from __future__ import division, print_function
import collections
import getopt
import h5py
import logging
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
import sys

matplotlib.rcParams['text.usetex'] = True
# Some of the labels used in the plotting routines require \mathscr
# Jump through some extra hoops to make sure we don't introduce
# duplicates into the preamble.
def unique_append(l, v):
    if v not in l:
        l.append(v)
unique_append(matplotlib.rcParams['text.latex.preamble'],
              r'\usepackage[charter]{mathdesign}')

# Prepare a logger to produce messages
logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
l = logging.getLogger(os.path.splitext(os.path.basename(__file__))[0])



class Bunch(dict):
    "An implementation of the Bunch pattern."
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self


def maybe_assume_uncorrelated(data, keyAB, keyA, keyB=None, warn=True):
    """
    If "AB" is not in dictionary data then assume Cov(A,B) = 0 implying
    E[AB] =  E[A] * E[B].  This is obviously not ideal from a rigor
    perspective and therefore warnings are emitted whenever warn is True.
    """
    if keyAB in data:
        pass
    elif keyB is None:
        data[keyAB] = data[keyA]**2
        l.warn("Estimated unknown %s by neglecting higher moment of %s"
               % (keyAB, keyA))
    else:
        data[keyAB] = data[keyA] * data[keyB]
        l.warn("Estimated unknown %s by assuming %s and %s are uncorrelated"
               % (keyAB, keyA, keyB))


class Data(object):
    """
    Storage of data loaded and computed from a single summary file.

    When enforceSymmetry is True and the file represents channel data,
    the upper and lower halves of the channel are assumed to be
    symmetric causing the profile data to be averaged.
    """


    # FIXME Correct symmetry enforcement for antisymmetric quantities
    def __init__(self, filename, enforceSymmetry=False):
        "Prepare averaged information and derived quantities from a filename"

        # Load file and prepare scalar and vector storage locations
        # Automatically unpack a variety of data mapping it into Bunches
        self.file    = h5py.File(filename, 'r')
        self.y       = self.file['y'][()]
        self.code    = Bunch(alpha = self.file['alpha'][0],
                             beta  = self.file['beta' ][0],
                             gamma = self.file['gamma'][0],
                             Ma    = self.file['Ma'   ][0],
                             Re    = self.file['Re'   ][0],
                             Pr    = self.file['Pr'   ][0])
        self.htdelta = self.file['htdelta'][0]

        # Additionally, "bl.*" or "chan.*" is mapped into self.*.
        self.bar    = Bunch()  # From "mu" attribute of bar_*
        self.bulk   = Bunch()  # Populated after construction
        self.rms    = Bunch()  # Populated after construction
        self.local  = Bunch()  # Populated after construction
        self.lower  = Bunch()  # From lower_*
        self.sigma  = Bunch()  # From "mu_sigma" attribute of bar_*
        self.tilde  = Bunch()  # Populated after construction
        self.tke    = Bunch()  # Populated after construction
        self.upper  = Bunch()  # From upper_*
        self.weight = Bunch()  # From *_weights
        for k, v in self.file.iteritems():
            if   k.startswith("bar_"):
                self.bar  [k[4:]] = v.attrs["mu"]
                self.sigma[k[4:]] = v.attrs["mu_sigma"]
            elif k.startswith("bulk_"):
                self.bulk[k[5:]] = v[0]
            elif k.startswith("lower_"):
                self.lower[k[6:]] = v[()]
            elif k.startswith("upper_"):
                self.upper[k[6:]] = v[()]
            elif k.endswith("_weights"):
                self.weight[k[:-8]] = v[()]
            elif k.startswith("bl."):
                self.__dict__[k[3:]] = Bunch(
                    **{a:b[0] for a,b in v.attrs.items()})
            elif k.startswith("chan."):
                self.__dict__[k[5:]] = Bunch(
                    **{a:b[0] for a,b in v.attrs.items()})

        # Employ symmetry assumptions, if requested and applicable
        # to improve mean estimates and associated uncertainties.
        # Based upon
        #     X ~ N(\mu_x, \sigma_x^2)
        #     Y ~ N(\mu_y, \sigma_y^2)
        #     Z = a X + b Y
        #       ~ N(a \mu_x + b \mu_y, a^2 \sigma_x^2 + b^2 \sigma_y^2)
        # implying
        #     std(Z) = sqrt(std(x)**2 + std(y)**2) / 2
        # which leads to the computations below.
        # Notice self.file is NOT adjusted, only mean profile information,
        # and also that the domain is not truncated to only the lower half.
        if self.htdelta >= 0 and enforceSymmetry:
            for k, v in self.bar.iteritems():
                self.bar[k] = (v + np.flipud(v)) / 2
            for k, v in self.sigma.iteritems():
                self.sigma[k] = np.sqrt(v**2 + np.flipud(v)**2) / 2
            self.enforceSymmetry = True
        else:
            self.enforceSymmetry = False

        # Compute scaling information per the incoming data, including...
        if "visc" in self.__dict__ and "wall" in self.__dict__:

            # ...classical incompressible plus units and...
            self.plus   = Bunch()
            self.plus.u = self.visc.u_tau
            self.plus.y = (self.wall.rho*self.y*self.plus.u / self.wall.mu
                          *self.code.Re)

            # ...variable density star units per Huang et al JFM 1995
            # (Beware this does weird things on the upper channel half)
            self.star   = Bunch()
            self.star.u = np.sqrt(self.visc.tau_w / self.bar.rho)
            self.star.y = (self.bar.rho*self.y*self.star.u / self.bar.mu
                          *self.code.Re)
        else:
            l.warn("Star and plus unit computations not performed"
                   " because /{bl,chan}.{visc,wall} not found")

        # Compute a whole mess of derived information, if possible
        try:
            from perfect_decl import pointwise
            pointwise(self.code.gamma,
                      self.code.Ma, self.code.Re, self.code.Pr, self.y,
                      self.bar, self.rms, self.tilde, self.local, self.tke)
        except ImportError:
            l.warn("Pointwise computations not performed"
                   " because module 'perfect_decl' not found")


    # FIXME Correct symmetry enforcement for antisymmetric quantities
    def esterrcov(self, *keys):
        """
        Produce an y-dependent empirical covariance matrix for the given
        keys returning squared standard errors from self.sigma on
        the diagonal.
        """
        # Check incoming arguments
        if not keys:
            raise ValueError("One or more keys must be provided")
        for k in keys:
            if k not in self.sigma:
                raise KeyError("Key %s not in sigma" % k)
            if "bar_"+k not in self.file:
                raise KeyError("Key %s not in sigma" % "bar_"+k)

        # Preallocate space to pack incoming data and store results
        ex = self.file["bar_"+k]
        nt, ny = ex.shape
        nv = len(keys)
        out = np.empty((nv, nv, ny), dtype=ex.dtype, order='C')
        dat = np.empty((nv, nt),     dtype=ex.dtype, order='C')

        # Process each y location in turn...
        for j in xrange(ny):
            # ...pack temporal trace into dat buffer
            # (if needed averaging with the opposite channel half)
            for i, key in enumerate(keys):
                dat[i, :] = self.file["bar_"+key][:,j]
            if self.enforceSymmetry:
                for i, key in enumerate(keys):
                    dat[i, :] += self.file["bar_"+key][:,ny-1-j]
                dat /= 2
            # ...obtain correlation matrix with unit diagonal
            t = np.corrcoef(dat, bias=1)
            # ...scale correlation by squared standard error
            # to obtain a covariance matrix.
            for i, key in enumerate(keys):
                t[i,:] *= self.sigma[key][j]
                t[:,i] *= self.sigma[key][j]
            # ...and store the result into out
            out[:, :, j] = t

        return out


def plot_profiles(data, fbottom=None, ftop=None, **fig_kw):
    """
    Plot mean primitive profiles, their RMS fluctuations, and uncertainties.
    """
    fig, ax = plt.subplots(2, 2,
                           sharex=('all' if (data.htdelta >= 0) else False),
                           squeeze=False, **fig_kw)
    bar, tilde, sigma, star = data.bar, data.tilde, data.sigma, data.star

    # Change our x axes slightly depending on channel versus boundary layer
    if data.htdelta >= 0:
        x  = star.y
        xl = None
    else:
        x  = data.y/data.thick.delta99
        xl = r"$y/\delta_{99}$"

    #########################################################################
    # Build dictionary of means and list of standard errors for upper row
    m = collections.OrderedDict()
    s = []

    m[r"$\bar{\rho}$"] = bar.rho
    s.append(sigma.rho)

    m[r"$\bar{u}$" ] = bar.u
    s.append(sigma.u)

    m[r"$\bar{T}$"] = bar.T
    s.append(sigma.T)

    m[r"$\bar{\mu}$"] = bar.mu
    s.append(sigma.mu)

    # Plot upper left subfigure
    for k, v in m.iteritems():
        ax[0][0].plot(x, v, label=k)
    if xl:
        ax[0][0].set_xlabel(xl)
    ax[0][0].set_ylabel(r"$\mu$")

    # Plot upper right subfigure
    for k, v in m.iteritems():
        ax[0][1].semilogy(x, s.pop(0)/v, label=k)
    if xl:
        ax[0][1].set_xlabel(xl)
    ax[0][1].set_ylabel(r"$\sigma_\mu / \mu$")
    if fbottom:
        ax[0][1].set_ylim(bottom=fbottom)
    if ftop:
        ax[0][1].set_ylim(top=ftop)

    #########################################################################
    # Build dictionary of means and standard errors for lower row
    # Variance expressions from postproc/propagation.py -d perfect.decl
    m.clear()
    del s[:]

    m[r"$\widetilde{u^{\prime\prime{}2}}$"] = tilde.upp_upp
    cov = data.esterrcov("rho",
                         "rho_u",
                         "rho_u_u")
    s.append(
        cov[0,0,:] * ((bar.rho*bar.rho_u_u - 2*bar.rho_u**2)**2/bar.rho**6)
      + cov[0,1,:] * (4*(bar.rho*bar.rho_u_u - 2*bar.rho_u**2)*bar.rho_u/bar.rho**5)
      + cov[0,2,:] * (2*(-bar.rho*bar.rho_u_u + 2*bar.rho_u**2)/bar.rho**4)
      + cov[1,1,:] * (4*bar.rho_u**2/bar.rho**4)
      + cov[1,2,:] * (-4*bar.rho_u/bar.rho**3)
      + cov[2,2,:] * (bar.rho**(-2))
    )

    m[r"$\widetilde{v^{\prime\prime{}2}}$"] = tilde.vpp_vpp
    cov = data.esterrcov("rho",
                         "rho_v",
                         "rho_v_v")
    s.append(
        cov[0,0,:] * ((bar.rho*bar.rho_v_v - 2*bar.rho_v**2)**2/bar.rho**6)
      + cov[0,1,:] * (4*(bar.rho*bar.rho_v_v - 2*bar.rho_v**2)*bar.rho_v/bar.rho**5)
      + cov[0,2,:] * (2*(-bar.rho*bar.rho_v_v + 2*bar.rho_v**2)/bar.rho**4)
      + cov[1,1,:] * (4*bar.rho_v**2/bar.rho**4)
      + cov[1,2,:] * (-4*bar.rho_v/bar.rho**3)
      + cov[2,2,:] * (bar.rho**(-2))
    )

    m[r"$\widetilde{w^{\prime\prime{}2}}$"] = tilde.wpp_wpp
    cov = data.esterrcov("rho",
                         "rho_w",
                         "rho_w_w")
    s.append(
        cov[0,0,:] * ((bar.rho*bar.rho_w_w - 2*bar.rho_w**2)**2/bar.rho**6)
      + cov[0,1,:] * (4*(bar.rho*bar.rho_w_w - 2*bar.rho_w**2)*bar.rho_w/bar.rho**5)
      + cov[0,2,:] * (2*(-bar.rho*bar.rho_w_w + 2*bar.rho_w**2)/bar.rho**4)
      + cov[1,1,:] * (4*bar.rho_w**2/bar.rho**4)
      + cov[1,2,:] * (-4*bar.rho_w/bar.rho**3)
      + cov[2,2,:] * (bar.rho**(-2))
    )

    m[r"$\widetilde{u^{\prime\prime}v^{\prime\prime}}$"] = tilde.upp_vpp
    cov = data.esterrcov("rho",
                         "rho_u",
                         "rho_u_v",
                         "rho_v")
    s.append(
        cov[0,0,:] * ((bar.rho*bar.rho_u_v - 2*bar.rho_u*bar.rho_v)**2/bar.rho**6)
      + cov[0,1,:] * (2*(bar.rho*bar.rho_u_v - 2*bar.rho_u*bar.rho_v)*bar.rho_v/bar.rho**5)
      + cov[0,2,:] * (2*(-bar.rho*bar.rho_u_v + 2*bar.rho_u*bar.rho_v)/bar.rho**4)
      + cov[0,3,:] * (2*(bar.rho*bar.rho_u_v - 2*bar.rho_u*bar.rho_v)*bar.rho_u/bar.rho**5)
      + cov[1,1,:] * (bar.rho_v**2/bar.rho**4)
      + cov[1,2,:] * (-2*bar.rho_v/bar.rho**3)
      + cov[1,3,:] * (2*bar.rho_u*bar.rho_v/bar.rho**4)
      + cov[2,2,:] * (bar.rho**(-2))
      + cov[2,3,:] * (-2*bar.rho_u/bar.rho**3)
      + cov[3,3,:] * (bar.rho_u**2/bar.rho**4)
    )

    m[r"$k$"] = tilde.k
    cov = data.esterrcov("rho",
                         "rho_u",
                         "rho_u_u",
                         "rho_v",
                         "rho_v_v",
                         "rho_w",
                         "rho_w_w")
    s.append(
      + cov[0,0,:] * ((bar.rho*bar.rho_u_u + bar.rho*bar.rho_v_v + bar.rho*bar.rho_w_w - 2*bar.rho_u**2 - 2*bar.rho_v**2 - 2*bar.rho_w**2)**2/(4*bar.rho**6))
      + cov[0,1,:] * ((bar.rho*bar.rho_u_u + bar.rho*bar.rho_v_v + bar.rho*bar.rho_w_w - 2*bar.rho_u**2 - 2*bar.rho_v**2 - 2*bar.rho_w**2)*bar.rho_u/bar.rho**5)
      + cov[0,2,:] * ((-bar.rho*bar.rho_u_u/2 - bar.rho*bar.rho_v_v/2 - bar.rho*bar.rho_w_w/2 + bar.rho_u**2 + bar.rho_v**2 + bar.rho_w**2)/bar.rho**4)
      + cov[0,3,:] * ((bar.rho*bar.rho_u_u + bar.rho*bar.rho_v_v + bar.rho*bar.rho_w_w - 2*bar.rho_u**2 - 2*bar.rho_v**2 - 2*bar.rho_w**2)*bar.rho_v/bar.rho**5)
      + cov[0,4,:] * ((-bar.rho*bar.rho_u_u/2 - bar.rho*bar.rho_v_v/2 - bar.rho*bar.rho_w_w/2 + bar.rho_u**2 + bar.rho_v**2 + bar.rho_w**2)/bar.rho**4)
      + cov[0,5,:] * ((bar.rho*bar.rho_u_u + bar.rho*bar.rho_v_v + bar.rho*bar.rho_w_w - 2*bar.rho_u**2 - 2*bar.rho_v**2 - 2*bar.rho_w**2)*bar.rho_w/bar.rho**5)
      + cov[0,6,:] * ((-bar.rho*bar.rho_u_u/2 - bar.rho*bar.rho_v_v/2 - bar.rho*bar.rho_w_w/2 + bar.rho_u**2 + bar.rho_v**2 + bar.rho_w**2)/bar.rho**4)
      + cov[1,1,:] * (bar.rho_u**2/bar.rho**4)
      + cov[1,2,:] * (-bar.rho_u/bar.rho**3)
      + cov[1,3,:] * (2*bar.rho_u*bar.rho_v/bar.rho**4)
      + cov[1,4,:] * (-bar.rho_u/bar.rho**3)
      + cov[1,5,:] * (2*bar.rho_u*bar.rho_w/bar.rho**4)
      + cov[1,6,:] * (-bar.rho_u/bar.rho**3)
      + cov[2,2,:] * (1/(4*bar.rho**2))
      + cov[2,3,:] * (-bar.rho_v/bar.rho**3)
      + cov[2,4,:] * (1/(2*bar.rho**2))
      + cov[2,5,:] * (-bar.rho_w/bar.rho**3)
      + cov[2,6,:] * (1/(2*bar.rho**2))
      + cov[3,3,:] * (bar.rho_v**2/bar.rho**4)
      + cov[3,4,:] * (-bar.rho_v/bar.rho**3)
      + cov[3,5,:] * (2*bar.rho_v*bar.rho_w/bar.rho**4)
      + cov[3,6,:] * (-bar.rho_v/bar.rho**3)
      + cov[4,4,:] * (1/(4*bar.rho**2))
      + cov[4,5,:] * (-bar.rho_w/bar.rho**3)
      + cov[4,6,:] * (1/(2*bar.rho**2))
      + cov[5,5,:] * (bar.rho_w**2/bar.rho**4)
      + cov[5,6,:] * (-bar.rho_w/bar.rho**3)
      + cov[6,6,:] * (1/(4*bar.rho**2))
    )

    # Plot lower left subfigure (includes normalization)
    for k, v in m.iteritems():
        ax[1][0].plot(star.y, v / star.u**2, label=k)
    ax[1][0].set_xlabel(r"$y^\ast$")
    ax[1][0].set_ylabel(r"$\mu / u_\tau^{\ast{}2}$")

    # Plot lower right subfigure (includes normalization)
    for k, v in m.iteritems():
        ax[1][1].semilogy(star.y, np.sqrt(s.pop(0))/np.abs(v), label=k)
    ax[1][1].set_xlabel(r"$y^\ast$")
    ax[1][1].set_ylabel(r"$\sigma_\mu / \left|\mu\right|$")
    if fbottom:
        ax[1][1].set_ylim(bottom=fbottom)
    if ftop:
        ax[1][1].set_ylim(top=ftop)

    if data.htdelta >= 0:
        # Truncate at half channel width
        ax[0][0].set_xlim(right=np.median(star.y))
        ax[0][1].set_xlim(right=np.median(star.y))
        ax[1][0].set_xlim(right=np.median(star.y))
        ax[1][1].set_xlim(right=np.median(star.y))
        # Add legends on rightmost images
        ax[0][1].legend(frameon=False, loc='best')
        ax[1][1].legend(frameon=False, loc='best')
    else:
        # Truncate at boundary layer thickness
        ax[0][0].set_xlim(right=1.0)
        ax[0][1].set_xlim(right=1.0)
        ax[1][0].set_xlim(right=star.y[np.argmax(data.y/data.thick.delta99 >= 1.0)])
        ax[1][1].set_xlim(right=star.y[np.argmax(data.y/data.thick.delta99 >= 1.0)])
        # Add legends on leftmost images
        ax[0][0].legend(frameon=False, loc='best')
        ax[1][0].legend(frameon=False, loc='best')

    return fig


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
# TODO Add a show_residual option like plot_rho, plot_rho_u, etc.
def plot_tke(data, y=None, vert=1, thresh=None, merge_pflux=False,
             ax=None, **plotargs):
    """
    Plot TKE budgets from data permitting rescaling and thresholding.

    If thresh is not None, it is taken as the fraction of maximum
    production used as a cutoff for suppressing lines.  Positive
    thresh retains only lines with maxima above that value while
    negative thresh retains only lines with maxima below that value.
    Guarini et al 2000 used 25.  Regardless, identically zero curves
    are always suppressed.

    If merge_pflux is True, these two pressure-related flux terms
         \bar{p}\nabla\cdot\overline{u''}/\mbox{Ma}^2$
        -\nabla\cdot\bar{\rho}\widetilde{T''u''}/\gamma/\mbox{Ma}^2$
    are reported as a single curve.
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Plot along y unless requested otherwise
    if y is None:
        y = data.y

    # Plotting cutoff based on magnitude of the production term.
    max_production = np.max(np.abs(vert * data.tke.production))
    def pthresh(y, q, *args, **kwargs):
        max_q = np.max(np.abs(q))
        if max_q == 0:
            pass # Always suppress identically zero data
        elif thresh is None:
            return ax.plot(y, q, *args, **kwargs)
        elif thresh >= 0:
            if max_q >= max_production / +thresh:
                return ax.plot(y, q, *args, **kwargs)
        elif thresh < 0:
            if max_q < max_production /  -thresh:
                return ax.plot(y, q, *args, **kwargs)
        else:
            assert False, "Sanity error"

    # Produce plots in order of most to least likely to exceed thresh
    # This causes any repeated linetypes to be fairly simple to distinguish
    # Likely to exceed threshold but we will check anyway
    pthresh(y, vert * data.tke.production,
            label=r"$- \bar{\rho}\widetilde{u''\otimes{}u''}:\nabla\tilde{u}$",
            **plotargs)
    pthresh(y, vert * data.tke.dissipation,
            label=r"$- \bar{\rho}\epsilon / \mbox{Re}$",
            **plotargs)
    pthresh(y, vert * data.tke.transport,
            label=r"$- \nabla\cdot \bar{\rho} \widetilde{{u''}^{2}u''} / 2$",
            **plotargs)
    pthresh(y, vert * data.tke.diffusion,
            label=r"$\nabla\cdot \overline{\tau{}u''}/\mbox{Re}$",
            **plotargs)
    # Conceivable that these will not exceed the threshold
    # Permit two different ways to view the pressure terms per merge_pflux
    if merge_pflux:
        pthresh(y,   vert * data.tke.pmassflux
                   + vert * data.tke.pheatflux,
                label=r"$\left("
                      r"\bar{p}\nabla\cdot\overline{u''}"
                      r"-\nabla\cdot\bar{\rho}\widetilde{T''u''}/\gamma"
                      r"\right)/\mbox{Ma}^2$",
                **plotargs)
        # Pressure dilatation appears inside conditional to preserve ordering
        pthresh(y, vert * data.tke.pdilatation,
                label=r"$\overline{p' \nabla\cdot{}u''}/\mbox{Ma}^2$",
                **plotargs)
    else:
        pthresh(y, vert * data.tke.pmassflux,
                label=r"$\bar{p}\nabla\cdot\overline{u''}/\mbox{Ma}^2$",
                **plotargs)
        # Pressure dilatation appears inside conditional to preserve ordering
        pthresh(y, vert * data.tke.pdilatation,
                label=r"$\overline{p' \nabla\cdot{}u''}/\mbox{Ma}^2$",
                **plotargs)
        pthresh(y, vert * data.tke.pheatflux,
                label=r"$-\nabla\cdot\bar{\rho}\widetilde{T''u''}/\gamma/\mbox{Ma}^2$",
                **plotargs)
    pthresh(y, vert * data.tke.convection,
            label=r"$- \nabla\cdot\bar{\rho}k\tilde{u}$",
            **plotargs)
    pthresh(y, vert * (data.tke.forcing + data.tke.constraint),
            label=r"$\overline{f\cdot{}u''}$",
            **plotargs)
    pthresh(y, vert * data.tke.slowgrowth,
            label=r"$\overline{\mathscr{S}_{\rho{}u}\cdot{}u''}$",
            **plotargs)

    return ax.figure


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
def plot_rho(data, y=None, vert=1, show_residual=False, ax=None, **plotargs):
    """
    Plot density budgets from data permitting rescaling.
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Plot along y unless requested otherwise
    if y is None:
        y = data.y

    # Gather all pointwise quantities into a label -> value dictionary
    curves = collections.OrderedDict()

    # See perfect.decl for more background
    curves[r"$- \nabla\cdot\bar{\rho}\tilde{u}$"] = (
       - data.bar.rho_v__y
    )
    curves[r"$  \overline{\mathscr{S}_{\rho}}$"] = (
       + data.bar.Srho
    )
    curves[r"$  \overline{\mathscr{C}_{\rho}}$"] = (
       + data.bar.Crho
    )

    # Produce the plot for nontrivial quantities
    if show_residual:
        residual = np.zeros_like(y)
    for key, val in curves.iteritems():
        if show_residual:
            residual += val
        if np.abs(val).max() > 0:
            ax.plot(y, vert * val, label=key, **plotargs)
    if show_residual:
        ax.plot(y, vert * residual, 'k:', label="Residual")

    return ax.figure


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
def plot_rho_u(data, y=None, vert=1, ax=None, show_residual=False, **plotargs):
    """
    Plot streamwise momentum budgets from data permitting rescaling.
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Plot along y unless requested otherwise
    if y is None:
        y = data.y

    # Gather all pointwise quantities into a label -> value dictionary
    curves = collections.OrderedDict()

    # See perfect.decl for more background
    curves[r"$- \nabla\cdot\left.\tilde{u}\otimes\bar{\rho}\tilde{u}\right.$"] = (
       - data.tilde.v*data.bar.rho_u__y - data.bar.rho_u*data.tilde.v__y
    )
    curves[r"$- \frac{1}{\textrm{Ma}^2}\nabla{}\bar{p}$"] = (
       - np.zeros_like(y) / data.code.Ma**2
    )
    curves[r"$  \nabla\cdot\left. \bar{\tau}/\textrm{Re} \right.$"] = (
       + data.bar.tauxy__y / data.code.Re
    )
    curves[r"$- \nabla\cdot\left. \bar{\rho} \widetilde{u''\otimes{}u''} \right.$"] = (
       - data.bar.rho*data.tilde.upp_vpp__y
       - data.tilde.upp_vpp*data.bar.rho__y
    )
    curves[r"$  \bar{f}$"] = (
       + data.bar.fx
    )
    curves[r"$  \overline{\mathscr{S}_{\rho{}u}}$"] = (
       + data.bar.Srhou
    )
    curves[r"$  \overline{\mathscr{C}_{\rho{}u}}$"] = (
       + data.bar.Crhou
    )

    # Produce the plot for nontrivial quantities
    if show_residual:
        residual = np.zeros_like(y)
    for key, val in curves.iteritems():
        if show_residual:
            residual += val
        if np.abs(val).max() > 0:
            ax.plot(y, vert * val, label=key, **plotargs)
    if show_residual:
        ax.plot(y, vert * residual, 'k:', label="Residual")

    return ax.figure


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
def plot_rho_v(data, y=None, vert=1, ax=None, show_residual=False, **plotargs):
    """
    Plot wall-normal momentum budgets from data permitting rescaling.
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Plot along y unless requested otherwise
    if y is None:
        y = data.y

    # Gather all pointwise quantities into a label -> value dictionary
    curves = collections.OrderedDict()

    # See perfect.decl for more background
    curves[r"$- \nabla\cdot\left.\tilde{u}\otimes\bar{\rho}\tilde{u}\right.$"] = (
       - data.tilde.v*data.bar.rho_v__y - data.bar.rho_v*data.tilde.v__y
    )
    curves[r"$- \frac{1}{\textrm{Ma}^2}\nabla{}\bar{p}$"] = (
       - data.bar.p__y / data.code.Ma**2
    )
    curves[r"$  \nabla\cdot\left. \bar{\tau}/\textrm{Re} \right.$"] = (
       + data.bar.tauyy__y / data.code.Re
    )
    curves[r"$- \nabla\cdot\left. \bar{\rho} \widetilde{u''\otimes{}u''} \right.$"] = (
       - data.bar.rho*data.tilde.vpp_vpp__y
       - data.tilde.vpp_vpp*data.bar.rho__y
    )
    curves[r"$  \bar{f}$"] = (
       + data.bar.fy
    )
    curves[r"$  \overline{\mathscr{S}_{\rho{}u}}$"] = (
       + data.bar.Srhov
    )
    curves[r"$  \overline{\mathscr{C}_{\rho{}u}}$"] = (
       + data.bar.Crhov
    )

    # Produce the plot for nontrivial quantities
    if show_residual:
        residual = np.zeros_like(y)
    for key, val in curves.iteritems():
        if show_residual:
            residual += val
        if np.abs(val).max() > 0:
            ax.plot(y, vert * val, label=key, **plotargs)
    if show_residual:
        ax.plot(y, vert * residual, 'k:', label="Residual")

    return ax.figure


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
def plot_rho_w(data, y=None, vert=1, ax=None, show_residual=False, **plotargs):
    """
    Plot spanwise momentum budgets from data permitting rescaling.
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Plot along y unless requested otherwise
    if y is None:
        y = data.y

    # Gather all pointwise quantities into a label -> value dictionary
    curves = collections.OrderedDict()

    # See perfect.decl for more background
    curves[r"$- \nabla\cdot\left.\tilde{u}\otimes\bar{\rho}\tilde{u}\right.$"] = (
       - data.tilde.v*data.bar.rho_w__y - data.bar.rho_w*data.tilde.v__y
    )
    curves[r"$- \frac{1}{\textrm{Ma}^2}\nabla{}\bar{p}$"] = (
       - np.zeros_like(y) / data.code.Ma**2
    )
    curves[r"$  \nabla\cdot\left. \bar{\tau}/\textrm{Re} \right.$"] = (
       + data.bar.tauyz__y / data.code.Re
    )
    curves[r"$- \nabla\cdot\left. \bar{\rho} \widetilde{u''\otimes{}u''} \right.$"] = (
       - data.bar.rho*data.tilde.vpp_wpp__y
       - data.tilde.vpp_wpp*data.bar.rho__y
    )
    curves[r"$  \bar{f}$"] = (
       + data.bar.fz
    )
    curves[r"$  \overline{\mathscr{S}_{\rho{}u}}$"] = (
       + data.bar.Srhow
    )
    curves[r"$  \overline{\mathscr{C}_{\rho{}u}}$"] = (
       + data.bar.Crhow
    )

    # Produce the plot for nontrivial quantities
    if show_residual:
        residual = np.zeros_like(y)
    for key, val in curves.iteritems():
        if show_residual:
            residual += val
        if np.abs(val).max() > 0:
            ax.plot(y, vert * val, label=key, **plotargs)
    if show_residual:
        ax.plot(y, vert * residual, 'k:', label="Residual")

    return ax.figure


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
def plot_rho_E(data, y=None, vert=1, ax=None, show_residual=False, **plotargs):
    """
    Plot total energy budgets from data permitting rescaling.
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Plot along y unless requested otherwise
    if y is None:
        y = data.y

    # Gather all pointwise quantities into a label -> value dictionary
    curves = collections.OrderedDict()

    # See perfect.decl for more background
    curves[r"$- \nabla\cdot\bar{\rho}\tilde{H}\tilde{u}$"] = (
       - data.bar.rho_v*data.tilde.H__y - data.tilde.H*data.bar.rho_v__y
    )
    curves[r"$  \textrm{Ma}^{2} \nabla\cdot\left. \frac{\bar{\tau} \tilde{u}}{\textrm{Re}} \right.$"] = (
       + (data.code.Ma**2/data.code.Re)*( data.tilde.u*data.bar.tauxy__y + data.bar.tauxy*data.tilde.u__y
                                        + data.tilde.v*data.bar.tauyy__y + data.bar.tauyy*data.tilde.v__y
                                        + data.tilde.w*data.bar.tauyz__y + data.bar.tauyz*data.tilde.w__y )
    )
    curves[r"$- \textrm{Ma}^{2} \nabla\cdot\left( \bar{\rho} \widetilde{u''\otimes{}u''} \right) \tilde{u}$"] = (
       - (data.code.Ma**2)*( data.bar.rho_u*data.tilde.upp_vpp__y+data.tilde.upp_vpp*data.bar.rho_u__y
                           + data.bar.rho_v*data.tilde.vpp_vpp__y+data.tilde.vpp_vpp*data.bar.rho_v__y
                           + data.bar.rho_w*data.tilde.vpp_wpp__y+data.tilde.vpp_wpp*data.bar.rho_w__y )
    )
    curves[r"$- \frac{1}{2}\textrm{Ma}^{2} \nabla\cdot\left. \bar{\rho}\widetilde{{u''}^{2}u''} \right.$"] = (
       - (data.code.Ma*data.code.Ma/2)*( data.bar.rho*data.tilde.upp2vpp__y + data.tilde.upp2vpp*data.bar.rho__y )
    )
    curves[r"$  \textrm{Ma}^{2} \nabla\cdot\left. \frac{\overline{\tau{}u''}}{\textrm{Re}} \right.$"] = (
       + (data.code.Ma*data.code.Ma/data.code.Re)*( data.bar.tauuppy__y )
    )
    curves[r"$  \frac{1}{\gamma-1} \nabla\cdot\left. \frac{ \bar{\mu} \widetilde{\nabla{}T} }{\textrm{Re}\textrm{Pr}} \right.$"] = (
       + (
            data.tilde.nu*data.bar.rho_grady_T__y + data.bar.rho_grady_T*data.tilde.nu__y
         ) / ((data.code.gamma-1)*data.code.Re*data.code.Pr)
    )
    curves[r"$  \frac{1}{\gamma-1} \nabla\cdot\left. \frac{ \bar{\rho} \widetilde{\nu'' \left(\nabla{}T\right)''} } { \textrm{Re}\textrm{Pr} } \right.$"] = (
       + (
            data.bar.rho*data.tilde.nupp_gradyTpp__y + data.tilde.nupp_gradyTpp*data.bar.rho__y
         ) / ((data.code.gamma-1)*data.code.Re*data.code.Pr)
    )
    curves[r"$- \frac{1}{\gamma-1} \nabla\cdot\left. \bar{\rho} \widetilde{T''u''} \right.$"] = (
       - (
            data.bar.rho*data.tilde.Tpp_vpp__y + data.tilde.Tpp_vpp*data.bar.rho__y
         ) / (data.code.gamma-1)
    )
    curves[r"$  \textrm{Ma}^{2} \bar{f}\cdot\tilde{u}$"] = (
       + data.code.Ma*data.code.Ma*(data.bar.fx*data.tilde.u + data.bar.fy*data.tilde.v + data.bar.fz*data.tilde.w)
    )
    curves[r"$  \textrm{Ma}^{2} \overline{f\cdot{}u''}$"] = (
       + data.code.Ma*data.code.Ma*data.bar.f_dot_upp
    )
    curves[r"$  \bar{q}_b$"] = (
       + data.bar.qb
    )
    curves[r"$  \overline{\mathscr{S}_{\rho{}E}}$"] = (
       + data.bar.SrhoE
    )
    curves[r"$  \overline{\mathscr{C}_{\rho{}E}}$"] = (
       + data.bar.CrhoE
    )

    # Produce the plot for nontrivial quantities
    if show_residual:
        residual = np.zeros_like(y)
    for key, val in curves.iteritems():
        if show_residual:
            residual += val
        if np.abs(val).max() > 0:
            ax.plot(y, vert * val, label=key, **plotargs)
    if show_residual:
        ax.plot(y, vert * residual, 'k:', label="Residual")

    return ax.figure


def traceframe(grepkey, fnames, zerotime=False):
    """
    Build a DataFrame from data within possibly many state.dat, qoi.dat,
    or bc.dat files using the 4th column 't' as the index.  Logging level,
    time since binary launch, and time step number information is omitted.
    If zerotime is True, time will be shifted so the first time
    index is zero.
    """
    r = pd.DataFrame()

    for fname in fnames:
        f = os.popen('grep "%s" %s' % (grepkey, fname))
        t = pd.read_table(f, index_col=3, sep=r"\s*")
        t = t.drop(t.columns[0:2]+t.columns[3:4], axis=1)
        r = r.combine_first(t)

    # Based on http://stackoverflow.com/questions/14110721
    # Inplace=True would be nice but it isn't available in older Pandas
    if zerotime:
        t = r.reset_index()
        t.t -= t.t[0]
        r = t.set_index(['t'])

    return r


# I wanted the newer matplotlib LogLocator so I stole it.  C'est dommage.
class LogLocator(ticker.Locator):
    """
    Determine the tick locations for log axes
    """

    def __init__(self, base=10.0, subs=[1.0], numdecs=4, numticks=15):
        """
        place ticks on the location= base**i*subs[j]
        """
        self.base(base)
        self.subs(subs)
        self.numticks = numticks
        self.numdecs = numdecs

    def base(self, base):
        """
        set the base of the log scaling (major tick every base**i, i integer)
        """
        self._base = base + 0.0

    def subs(self, subs):
        """
        set the minor ticks the log scaling every base**i*subs[j]
        """
        if subs is None:
            self._subs = None  # autosub
        else:
            self._subs = np.asarray(subs) + 0.0

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        b = self._base
        # dummy axis has no axes attribute
        if hasattr(self.axis, 'axes') and self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax) / math.log(b))
            decades = np.arange(vmax - self.numdecs, vmax)
            ticklocs = b ** decades

            return ticklocs

        if vmin <= 0.0:
            if self.axis is not None:
                vmin = self.axis.get_minpos()

            if vmin <= 0.0 or not np.isfinite(vmin):
                raise ValueError(
                    "Data has no positive values, and therefore can not be "
                    "log-scaled.")

        vmin = math.log(vmin) / math.log(b)
        vmax = math.log(vmax) / math.log(b)

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        numdec = math.floor(vmax) - math.ceil(vmin)

        if self._subs is None:  # autosub
            if numdec > 10:
                subs = np.array([1.0])
            elif numdec > 6:
                subs = np.arange(2.0, b, 2.0)
            else:
                subs = np.arange(2.0, b)
        else:
            subs = self._subs

        stride = 1
        while numdec / stride + 1 > self.numticks:
            stride += 1

        decades = np.arange(math.floor(vmin) - stride,
                            math.ceil(vmax) + 2 * stride, stride)
        if hasattr(self, '_transform'):
            ticklocs = self._transform.inverted().transform(decades)
            if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
                ticklocs = np.ravel(np.outer(subs, ticklocs))
        else:
            if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
                ticklocs = []
                for decadeStart in b ** decades:
                    ticklocs.extend(subs * decadeStart)
            else:
                ticklocs = b ** decades

        return self.raise_if_exceeds(np.asarray(ticklocs))

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'
        b = self._base

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax) / math.log(b))
            vmin = b ** (vmax - self.numdecs)
            return vmin, vmax

        minpos = self.axis.get_minpos()

        if minpos <= 0 or not np.isfinite(minpos):
            raise ValueError(
                "Data has no positive values, and therefore can not be "
                "log-scaled.")

        if vmin <= minpos:
            vmin = minpos

        if not is_decade(vmin, self._base):
            vmin = decade_down(vmin, self._base)
        if not is_decade(vmax, self._base):
            vmax = decade_up(vmax, self._base)

        if vmin == vmax:
            vmin = decade_down(vmin, self._base)
            vmax = decade_up(vmax, self._base)
        result = mtransforms.nonsingular(vmin, vmax)
        return result


# TODO Split into two plots-- scenario tracking vs relaminarization
# TODO Add Re_\tau as a precursor to fluctuation collapse
def plot_relaminarization(dnames,
                          Re_theta=None,
                          Ma_e=None,
                          ratio_T=None,
                          v_wallplus=None,
                          p_ex=None,
                          delta99=None,
                          zerotime=False,
                          **kwargs):
    """
    Prepare a relaminarization study plot given job directories dnames.
    If provided, Ma_e, p_ex, etc. are used to plot target values.
    """
    # Load the data from various source files
    pbulk = traceframe('prod.bulk',  (s+"/qoi.dat" for s in dnames), zerotime)
    pg    = traceframe('bl.pg',      (s+"/qoi.dat" for s in dnames), zerotime)
    qoi   = traceframe('bl.qoi',     (s+"/qoi.dat" for s in dnames), zerotime)
    Re    = traceframe('bl.Re',      (s+"/qoi.dat" for s in dnames), zerotime)
    thick = traceframe('bl.thick',   (s+"/qoi.dat" for s in dnames), zerotime)
    visc  = traceframe('bl.visc',    (s+"/qoi.dat" for s in dnames), zerotime)
    famax = traceframe('favre.amax', (s+"/qoi.dat" for s in dnames), zerotime)

    # Produce a relaminarization summary figure
    fig = plt.figure(**kwargs)
    #
    ax = fig.add_subplot(911)
    ax.ticklabel_format(useOffset=False)
    ax.plot(Re.index, Re['Re_delta2'].values)
    if Re_theta:
        ax.hlines(Re_theta, qoi.index.min(), qoi.index.max(), 'r', '-.')
    ax.set_ylabel(r'$\mbox{Re}_{\theta}$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(912)
    ax.ticklabel_format(useOffset=False)
    ax.plot(qoi.index, qoi['Ma_e'].values)
    if Ma_e:
        ax.hlines(Ma_e, qoi.index.min(), qoi.index.max(), 'r', '-.')
    ax.set_ylabel(r'$\mbox{Ma}_{e}$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(913)
    ax.ticklabel_format(useOffset=False)
    ax.plot(qoi.index, qoi['ratio_T'].values)
    if ratio_T:
        ax.hlines(ratio_T, qoi.index.min(), qoi.index.max(), 'r', '-.')
    ax.set_ylabel(r'$T_e/T_w$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(914)
    ax.ticklabel_format(useOffset=False)
    ax.plot(visc.index, visc['v_wallplus'].values)
    if v_wallplus:
        ax.hlines(v_wallplus, visc.index.min(), visc.index.max(), 'r', '-.')
    ax.set_ylabel(r'$v_w^+$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(915)
    ax.ticklabel_format(useOffset=False)
    ax.plot(pg.index, pg['p_ex'].values)
    if p_ex:
        ax.hlines(p_ex, pg.index.min(), pg.index.max(), 'r', '-.')
    ax.set_ylabel(r'$p^\ast_{e,\xi}$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(916)
    ax.ticklabel_format(useOffset=False)
    ax.plot(thick.index, thick['delta99'].values)
    if delta99:
        ax.hlines(delta99, thick.index.min(), thick.index.max(), 'r', '-.')
    ax.set_ylabel(r'$\delta_{99}$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(917)
    ax.ticklabel_format(useOffset=False)
    ax.plot(pbulk.index, pbulk['total'].values)
    ax.set_ylabel("Integrated\n"
                  #r"$\overline{\rho u'' \otimes{} u''}:\nabla\tilde{u}$",
                  "production",
                  multialignment='center')
    if np.log10(pbulk['total'].values.max() / pbulk['total'].values.min()) > 1:
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(LogLocator(numticks=3))
    else:
        ax.set_yscale('linear')
        ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    if ax.get_ylim()[0] < 1e-9:
        ax.set_ylim(bottom=1e-9)
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(918)
    ax.ticklabel_format(useOffset=False)
    ax.plot(famax.index, np.abs(famax['uu'].values))
    ax.plot(famax.index, np.abs(famax['uv'].values))
    ax.plot(famax.index, np.abs(famax['uw'].values))
    ax.plot(famax.index, np.abs(famax['vv'].values))
    ax.plot(famax.index, np.abs(famax['vw'].values))
    ax.plot(famax.index, np.abs(famax['ww'].values))
    ax.set_yscale('log')
    ax.set_ylabel("max\n"
                  r"$\left|\widetilde{u_i''u_j''}\right|$",
                  multialignment='center')
    if ax.get_ylim()[0] < 1e-8:
        ax.set_ylim(bottom=1e-8)
    ax.yaxis.set_major_locator(LogLocator(numticks=3))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(919)
    ax.ticklabel_format(useOffset=False)
    ax.plot(visc.index, visc['cf'].values)
    ax.set_ylabel(r'$c_f$')
    ax.yaxis.set_major_locator(ticker.MaxNLocator( 4))
    #
    ax.set_xlabel(r'Nondimensional time $t\,u_0 / l_0$')
    #
    for ax in fig.axes:
        ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
        ax.margins(0.025, 0.08)
    #
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.10)

    return fig


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    outsuffix = None
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print(__doc__)
                return 0
            elif o == "-o":
                outsuffix = a
    except Usage as err:
        print(err.msg, file=sys.stderr)
        return 2

    # Push interactive mode off (in case we get used from IPython)
    was_interactive = plt.isinteractive()
    plt.interactive(False)

    # If not saving, then display.
    if not outsuffix:
        plt.show()

    # Pop interactive mode
    plt.interactive(was_interactive)


if __name__ == "__main__":
    sys.exit(main())
