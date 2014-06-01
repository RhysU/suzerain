#!/usr/bin/env python
"""Usage: summary_surf.py    [OPTIONS...] H5SUMMARY DATASET...
          summary_surf.py -f [OPTIONS...] H5SUMMARY DATASET1 DATASET2 DATASET3
The first form produces surfaces for each named DATASET in the file H5SUMMARY.
The second form plots DATASET1 - DATASET2*DATASET3 which can show fluctuations.

Options:
    -c CSTRIDE    Downsample by providing cstride=CSTRIDE to plot_surface
    -f            Plot fluctuating quantities using the second invocation type
    -h            Display this help message and exit
    -l LINEWIDTH  Set a non-zero linewidth=LINEWIDTH to plot_surface
    -o OUTSUFFIX  Save the output file DATASET.OUTSUFFIX instead of displaying
    -r RSTRIDE    Downsample by providing rstride=RSTRIDE to plot_surface
    -t TITLE      Set the specified plot title
    -C LEVELS     Instead of drawing surfaces, plot LEVELS contours.
    -T LO,HI      Limit plot extents to /t in [LO, HI]
    -Y LO,HI      Limit plot extents to /y in [LO, HI]
    -Z LO,HI      Limit plot extents to data values in [LO, HI]

File H5SUMMARY must have 1D /y and /t datasets to provide the meshgrid.
Each DATASET must have extents matching the /y and /t dataset sizes.
"""
from __future__ import print_function
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


def surface(x, y, xyplane, kwargs={}):
    fig  = plt.figure()
    ax   = fig.gca(projection='3d')
    X, Y = np.meshgrid(x, y)
    plot = ax.plot_surface(X, Y, xyplane,
                           vmin=np.min(xyplane), vmax=np.max(xyplane),
                           **kwargs)
    if 'cmap' in kwargs:
        cbar = fig.colorbar(plot)
    else:
        cbar = None
    return (fig, ax, plot, cbar)


def contour(x, y, xyplane, levels, kwargs={}):
    fig  = plt.figure()
    ax   = fig.gca()
    X, Y = np.meshgrid(x, y)
    plot = ax.contourf(X, Y, xyplane, levels, **kwargs)
    if 'cmap' in kwargs:
        cbar = fig.colorbar(plot)
    else:
        cbar = None
    return (fig, ax, plot, cbar)


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    cstride   = 1
    fluct     = False
    linewidth = 0
    outsuffix = None
    rstride   = 1
    levels    = 0
    title     = None
    textents  = (-np.inf, np.inf)
    yextents  = (-np.inf, np.inf)
    zextents  = None
    try:
        try:
            opts, args = getopt.getopt(argv[1:],
                                       "c:fhl:o:r:t:C:T:Y:Z:", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if   o == "-c":
                cstride = int(a)
            if   o == "-f":
                fluct = True
            elif o in ("-h", "--help"):
                print(__doc__)
                return 0
            elif o == "-l":
                linewidth = int(a)
            elif o == "-o":
                outsuffix = a
            elif o == "-r":
                rstride = int(a)
            elif o == "-t":
                title = a
            elif o == "-C":
                levels = int(a)
            elif o == "-T":
                textents = tuple(float(r) for r in a.split(','))
            elif o == "-Y":
                yextents = tuple(float(r) for r in a.split(','))
            elif o == "-Z":
                zextents = tuple(float(r) for r in a.split(','))
        if fluct:
            if len(args) != 4:
                print("Incorrect argument count.  See --help.", file=sys.stderr)
                return 1
        else:
            if len(args) < 2:
                print("Too few arguments.  See --help.", file=sys.stderr)
                return 1
    except Usage as err:
        print(err.msg, file=sys.stderr)
        return 2

    # Push interactive mode off (in case we get used from IPython)
    was_interactive = plt.isinteractive()
    plt.interactive(False)

    # Open the HDF5 file and read extent information
    # (Attempting to trim loaded data to be as small as possible)
    h5file = h5py.File(args[0], 'r')
    y      = h5file['/y'][()]
    yb     = np.nonzero(y > yextents[0])[0][ 0]
    ye     = np.nonzero(y < yextents[1])[0][-1]
    y      = y[yb:ye+1]
    t      = h5file['/t'][()]
    tb     = np.nonzero(t > textents[0])[0][ 0]
    te     = np.nonzero(t < textents[1])[0][-1]
    t      = t[tb:te+1]

    # Load everything, truncating any irrelevant range or values
    data, name = [], []
    for arg in args[1:]:
        data.append(h5file[arg][tb:te+1, yb:ye+1])
        name.append(arg)
        if zextents:
            data[-1][data[-1] < zextents[0]] = np.nan
            data[-1][data[-1] > zextents[1]] = np.nan

    # Compute derived fluctuating dataset, if requested
    if fluct:
        data3, name3 = data.pop(), name.pop()
        data2, name2 = data.pop(), name.pop()
        data1, name1 = data.pop(), name.pop()
        assert len(data) == 0
        assert len(name) == 0
        data.append(data1 - data2*data3)
        name.append("%s-%s*%s" % (name1, name2, name3))

    # Plot each dataset, annotating, and possibly saving
    for dataname, dataset in zip(name, data):

        if (levels == 0):
            fig, ax, plot, cbar = surface(y, t, dataset, {
                                          'cstride'  : cstride
                                        , 'rstride'  : rstride
                                        , 'cmap'     : matplotlib.cm.RdYlBu_r
                                        , 'linewidth': linewidth
                                        })
            ax.set_xlabel('Wall-normal distance')
            ax.set_ylabel('Simulation time')
            ax.set_zlabel(dataname)
        else:
            fig, ax, plot, cbar = contour(t, y, dataset.transpose(), levels, {
                                          'cmap'     : matplotlib.cm.RdYlBu_r
                                        , 'linewidth': linewidth
                                        })
            ax.set_xlabel('Simulation time')
            ax.set_ylabel('Wall-normal distance')
        if title:
            ax.set_title(title)
        if outsuffix:
            fig.savefig(dataname+'.'+outsuffix)
            plt.close(fig)

    # If not saving, then display.
    if not outsuffix:
        plt.show()

    # Pop interactive mode
    plt.interactive(was_interactive)


if __name__ == "__main__":
    sys.exit(main())
