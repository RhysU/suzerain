#!/usr/bin/env python
"""Usage: summary_surf.py [OPTIONS...] H5SUMMARY DATASET...
Produce surface plots for each named DATASET in the file H5SUMMARY.

Options:
    -c CSTRIDE    Downsample by providing cstride=CSTRIDE to plot_surface
    -h            Display this help message and exit
    -l LINEWIDTH  Set a non-zero linewidth=LINEWIDTH to plot_surface
    -o OUTSUFFIX  Save the output file DATASET.OUTSUFFIX instead of displaying
    -r RSTRIDE    Downsample by providing rstride=RSTRIDE to plot_surface
    -t TITLE      Set the specified plot title

File H5SUMMARY must have 1D /y and /t datasets to provide the meshgrid.
Each DATASET must have extents matching the /y and /t dataset sizes.
"""
import getopt
import h5py
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import numpy as np
import sys

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def surface(y, t, ytplane, kwargs={}):
    fig = plt.figure()
    ax  = fig.gca(projection='3d')
    Y,T = np.meshgrid(y, t)
    surf = ax.plot_surface(Y, T, ytplane,
                           vmin=np.min(ytplane), vmax=np.max(ytplane),
                           **kwargs)
    if 'cmap' in kwargs:
        cbar = fig.colorbar(surf)
    else:
        cbar = None
    return (fig, ax, surf, cbar)

def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    cstride   = 1
    linewidth = 0
    outsuffix = None
    rstride   = 1
    title     = None
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "c:hl:o:r:t:", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if   o == "-c":
                cstride = int(a)
            elif o in ("-h", "--help"):
                print __doc__
                return 0
            elif o == "-l":
                linewidth = int(a)
            elif o == "-o":
                outsuffix = a
            elif o == "-r":
                rstride = int(a)
            elif o == "-t":
                title = a
        if len(args) < 2:
            print >>sys.stderr, "Too few arguments.  See --help."
            return 1
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Push interactive mode off (in case we get used from IPython)
    was_interactive = plt.isinteractive()
    plt.interactive(False)

    # Open the HDF5 file and plot each dataset in turn
    h5file = h5py.File(args[0], 'r')
    y      = h5file['/y']
    t      = h5file['/t']
    for dataset in args[1:]:
        fig, ax, surf, cbar = surface(y, t, h5file[dataset], {
                                       'cstride'  : cstride
                                     , 'rstride'  : rstride
                                     , 'cmap'     : matplotlib.cm.RdYlBu_r
                                     , 'linewidth': linewidth
                                     })
        ax.set_xlabel('Wall-normal distance')
        ax.set_ylabel('Simulation time')
        if title:
            ax.set_title(title)
        if outsuffix:
            fig.savefig(dataset+'.'+outsuffix)
            plt.close(fig)

    if not outsuffix:
        plt.show()

    # Pop interactive mode
    plt.interactive(was_interactive)

if __name__ == "__main__":
    sys.exit(main())
