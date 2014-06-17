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


def evolution(pat, fnames):
    """
    Build a DataFrame from data within possibly many state.dat, qoi.dat, or
    bc.dat files using the 4th column 't' as the index.  Logging level,
    time since binary launch, and time step number information is omitted.
    """
    r = pd.DataFrame()
    for fname in fnames:
        f = os.popen('grep "%s" %s' % (pat, fname))
        t = pd.read_table(f, index_col=3, sep=r"\s*")
        t = t.drop(t.columns[0:2]+t.columns[3:4], axis=1)
        r = r.combine_first(t)
    return r


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
