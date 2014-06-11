#!/usr/bin/env python
"""Usage: plot_perfect.py [OPTIONS...] H5FILE...
Produces TBD

Options:
    -h            Display this help message and exit
    -o OUTSUFFIX  Save the output file PLOT.OUTSUFFIX instead of displaying

File H5FILE should have been produced by the perfect_summary application.
"""
from __future__ import division, print_function
import getopt
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys


class Bunch(dict):
    "An implementation of the Bunch pattern."
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self


class Mean(object):
    "Storage of data loaded and computed from a single summary file."
    def __init__(self, filename):
        "Load relevant Reynolds averaged information from filename"

        # Load file and prepare scalar and vector storage locations
        self.file  = h5py.File(filename, 'r')
        self.alpha = self.file['alpha'][0]
        self.beta  = self.file['beta' ][0]
        self.gamma = self.file['gamma'][0]
        self.Ma    = self.file['Ma'   ][0]
        self.Re    = self.file['Re'   ][0]
        self.Pr    = self.file['Pr'   ][0]
        self.y     = self.file['y'    ][()]

        # Automatically unpack a variety of data mapping it into Bunches
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
