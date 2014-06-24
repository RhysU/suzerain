#!/usr/bin/env python
"""Usage: perfect.py [OPTIONS...] H5FILE...
Produces TBD

Options:
    -h            Display this help message and exit
    -o OUTSUFFIX  Save the output file PLOT.OUTSUFFIX instead of displaying

File H5FILE should have been produced by the perfect_summary application.
"""
from __future__ import division, print_function
import getopt
import h5py
import logging
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys


# Prepare a logger to produce progress messages
logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
l = logging.getLogger(os.path.splitext(os.path.basename(__file__))[0])


class Bunch(dict):
    "An implementation of the Bunch pattern."
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self


class Data(object):
    "Storage of data loaded and computed from a single summary file."
    def __init__(self, filename):
        "Prepare averaged information and derived quantities from a filename"

        # Load file and prepare scalar and vector storage locations
        # Automatically unpack a variety of data mapping it into Bunches
        self.file = h5py.File(filename, 'r')
        self.y    = self.file['y'][()]
        self.code = Bunch(alpha = self.file['alpha'][0],
                          beta  = self.file['beta' ][0],
                          gamma = self.file['gamma'][0],
                          Ma    = self.file['Ma'   ][0],
                          Re    = self.file['Re'   ][0],
                          Pr    = self.file['Pr'   ][0])

        # Additionally, "bl.*" or "chan.*" is mapped into self.*.
        self.bar    = Bunch()  # From "mu" attribute of bar_*
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

        # Compute scaling information per the incoming data, including...
        if "visc" in self.__dict__ and "wall" in self.__dict__:

            # ...classical incompressible plus units and...
            self.plus   = Bunch()
            self.plus.u = self.visc.u_tau
            self.plus.y = (self.wall.rho*self.y*self.plus.u / self.wall.mu
                          *self.code.Re)

            # ...variable density star units per Huang et al JFM 1995
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
                      self.bar, self.tilde, self.local, self.tke)
        except ImportError:
            l.warn("Pointwise computations not performed"
                   " because module 'perfect_decl' not found")


# TODO Smooth per B-splines using ' from scipy.interpolate import interp1d'
def plot_tke(data, horz=1, vert=1, thresh=25, ax=None, **plotargs):
    """
    Plot TKE budgets from data permitting rescaling and thresholding
    """

    # Get a new axis if one was not supplied
    if not ax:
        fig, ax = plt.subplots()

    # Rescale the data as requested, permitting plus or star units.
    # Merging constraint and forcing as they're identical physically
    # and their separation is purely an artifact of the implementation.
    y           = horz * data.y
    convection  = vert * data.tke.convection
    production  = vert * data.tke.production
    dissipation = vert * data.tke.dissipation
    transport   = vert * data.tke.transport
    diffusion   = vert * data.tke.diffusion
    pmassflux   = vert * data.tke.pmassflux
    pdilatation = vert * data.tke.pdilatation
    pheatflux   = vert * data.tke.pheatflux
    forcing     = vert * (data.tke.forcing + data.tke.constraint)
    slowgrowth  = vert * data.tke.slowgrowth

    # Plotting cutoff based on magnitude of the production term
    # following Guarini et al JFM 2000 page 23.
    thresh = np.max(np.abs(production)) / thresh

    # Produce plots in order of most to least likely to exceed thresh
    # This causes any repeated linetypes to be fairly simple to distinguish
    def pthresh(y, q, *args, **kwargs):
        if np.max(np.abs(q)) > thresh:
            ax.plot(y, q, *args, **kwargs)
    # Likely to exceed threshold but we will check anyway
    pthresh(y, production, linestyle='-',
            label=r"$- \bar{\rho}\widetilde{u''\otimes{}u''}:\nabla\tilde{u}$",
            **plotargs)
    pthresh(y, dissipation, linestyle='-',
            label=r"$- \bar{\rho}\epsilon / \mbox{Re}$",
            **plotargs)
    pthresh(y, transport, linestyle='--',
            label=r"$- \nabla\cdot \bar{\rho} \widetilde{{u''}^{2}u''} / 2$",
            **plotargs)
    pthresh(y, diffusion, linestyle='-.',
            label=r"$\nabla\cdot \overline{\tau{}u''}/\mbox{Re}$",
            **plotargs)
    pthresh(y, slowgrowth, linestyle=':',
            label=r"$\overline{\mathscr{S}_{\rho{}u}\cdot{}u''}$",
            **plotargs)
    # Conceivable that these will not exceed the threshold
    pthresh(y, pmassflux, linestyle='-',
            label=r"$\bar{p}\nabla\cdot\overline{u''}/\mbox{Ma}^2$",
            **plotargs)
    pthresh(y, pdilatation, linestyle='--',
            label=r"$\overline{p' \nabla\cdot{}u''}/\mbox{Ma}^2$",
            **plotargs)
    pthresh(y, pheatflux, linestyle='-.',
            label=r"$-\nabla\cdot\bar{\rho}\widetilde{T''u''}/\gamma/\mbox{Ma}^2$",
            **plotargs)
    pthresh(y, convection, linestyle=':',
            label=r"$- \nabla\cdot\bar{\rho}k\tilde{u}$",
            **plotargs)
    pthresh(y, forcing, linestyle=':',
            label=r"$\overline{f\cdot{}u''}$",
            **plotargs)

    return ax


def traceframe(grepkey, fnames):
    """
    Build a DataFrame from data within possibly many state.dat, qoi.dat, or
    bc.dat files using the 4th column 't' as the index.  Logging level,
    time since binary launch, and time step number information is omitted.
    """
    r = pd.DataFrame()
    for fname in fnames:
        f = os.popen('grep "%s" %s' % (grepkey, fname))
        t = pd.read_table(f, index_col=3, sep=r"\s*")
        t = t.drop(t.columns[0:2]+t.columns[3:4], axis=1)
        r = r.combine_first(t)
    return r


# TODO Split into two plots-- scenario tracking vs relaminarization
# TODO Add Re_\tau as a precursor to fluctuation collapse
def plot_relaminarization(dnames,
                          Re_theta=None,
                          Ma_e=None,
                          ratio_T=None,
                          v_wallplus=None,
                          p_ex=None,
                          delta99=None,
                          **kwargs):
    """
    Prepare a relaminarization study plot given job directories dnames.
    If provided, Ma_e, p_ex, etc. are used to plot target values.
    """
    # Load the data from various source files
    pbulk = traceframe('prod.bulk',  (s+"/qoi.dat" for s in dnames))
    pg    = traceframe('bl.pg',      (s+"/qoi.dat" for s in dnames))
    qoi   = traceframe('bl.qoi',     (s+"/qoi.dat" for s in dnames))
    Re    = traceframe('bl.Re',      (s+"/qoi.dat" for s in dnames))
    thick = traceframe('bl.thick',   (s+"/qoi.dat" for s in dnames))
    visc  = traceframe('bl.visc',    (s+"/qoi.dat" for s in dnames))
    famax = traceframe('favre.amax', (s+"/qoi.dat" for s in dnames))

    # Produce a relaminarization summary figure
    fig = plt.figure(**kwargs)
    def yinclude(axis, val):
        ymin, ymax = axis.get_ylim()
        ymin = min(ymin, val - (ymax - ymin)/50)
        ymax = max(ymax, val + (ymax - ymin)/50)
        axis.set_ylim([ymin, ymax])
    #
    ax = fig.add_subplot(911)
    ax.ticklabel_format(useOffset=False)
    ax.plot(Re.index, Re['Re_delta2'].values)
    if Re_theta:
        yinclude(ax, Re_theta)
        ax.hlines(Re_theta, qoi.index.min(), qoi.index.max(), 'r', '-.')
    ax.set_ylabel(r'$\mbox{Re}_{\theta}$')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(912)
    ax.ticklabel_format(useOffset=False)
    ax.plot(qoi.index, qoi['Ma_e'].values)
    if Ma_e:
        yinclude(ax, Ma_e)
        ax.hlines(Ma_e, qoi.index.min(), qoi.index.max(), 'r', '-.')
    ax.set_ylabel(r'$\mbox{Ma}_{e}$')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(913)
    ax.ticklabel_format(useOffset=False)
    ax.plot(qoi.index, qoi['ratio_T'].values)
    if ratio_T:
        yinclude(ax, ratio_T)
        ax.hlines(ratio_T, qoi.index.min(), qoi.index.max(), 'r', '-.')
    ax.set_ylabel(r'$T_e/T_w$')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(914)
    ax.ticklabel_format(useOffset=False)
    ax.plot(visc.index, visc['v_wallplus'].values)
    if v_wallplus:
        yinclude(ax, v_wallplus)
        ax.hlines(v_wallplus, visc.index.min(), visc.index.max(), 'r', '-.')
    ax.set_ylabel(r'$v_w^+$')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(915)
    ax.ticklabel_format(useOffset=False)
    ax.plot(pg.index, pg['p_ex'].values)
    if p_ex:
        yinclude(ax, p_ex)
        ax.hlines(p_ex, pg.index.min(), pg.index.max(), 'r', '-.')
    ax.set_ylabel(r'$p^\ast_{e,\xi}$')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(916)
    ax.ticklabel_format(useOffset=False)
    ax.plot(thick.index, thick['delta99'].values)
    if delta99:
        yinclude(ax, delta99)
        ax.hlines(delta99, thick.index.min(), thick.index.max(), 'r', '-.')
    ax.set_ylabel(r'$\delta_{99}$')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(917)
    ax.ticklabel_format(useOffset=False)
    ax.plot(pbulk.index, pbulk['total'].values)
    ax.set_yscale('log')
    ax.set_ylabel("bulk\n"
                  #r"$\overline{\rho u'' \otimes{} u''}:\nabla\tilde{u}$",
                  "production",
                  multialignment='center')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
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
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    #
    ax = fig.add_subplot(919)
    ax.ticklabel_format(useOffset=False)
    ax.plot(visc.index, visc['cf'].values)
    ax.set_ylabel(r'cf')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator( 4))
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(10))
    #
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.20)

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
