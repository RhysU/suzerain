#!/usr/bin/env python
"""Usage: plot_stats.py HDF5FILE
Plot wall-normal profiles averaged over (X,Z) directions from HDF5FILE.
"""

# TODO: Add ability to plot things other than bar_rho.

# TODO: Add ability to take multiple h5files and plot all on same
# axes.  Note that the main supports this but plot fcn just
# overwrites.

import sys
import getopt
import h5py
import numpy as np
import gb 
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def plot(hdf5file):
    print "Plotting", hdf5file

    # Load a stats file
    f = h5py.File(hdf5file,'r')
        
    # Grab number of collocation points and B-spline order
    Ny=f['Ny'].value
    k=f['k'].value
    
    # Grab collocation points
    y = f['collocation_points_y'].value

    # Get "mass" matrix and convert to dense format
    D0T_gb = f['Dy0T'].value
    D0T = gb.gb2ge(D0T_gb, Ny, k-2)

    # Grab rho coefficients
    rho_coeff = f['bar_rho'].value
    rho_coeff = np.array(rho_coeff).reshape(Ny,1)


    # Grab rho_u coefficients
    rho_u_coeff = f['bar_rho_u'].value
    rho_u_coeff = np.array(rho_u_coeff).transpose().reshape(Ny,3)

    # Grab rho_E coefficients
    rho_E_coeff = f['bar_rho_E'].value
    rho_E_coeff = np.array(rho_E_coeff).reshape(Ny,1)

    # Done getting data
    f.close()
    
    D0 = D0T.transpose()

    # Coefficients -> Collocation points
    rho_col   = D0*rho_coeff
    rho_u_col = D0*rho_u_coeff
    rho_E_col = D0*rho_E_coeff


    print rho_u_col
    print rho_u_col[:,0]

    # Plots
    pyplot.figure(1)
    rho_key = "bar_rho_" + str(0)
    pyplot.plot(y, rho_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho.eps', bbox_inches='tight')

    pyplot.figure(2)
    rho_key = "bar_rho_u" + str(0)
    pyplot.plot(y, rho_u_col[:,0], linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_u.eps', bbox_inches='tight')

    pyplot.figure(3)
    rho_key = "bar_rho_v" + str(0)
    pyplot.plot(y, rho_u_col[:,1], linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_v.eps', bbox_inches='tight')

    pyplot.figure(4)
    rho_key = "bar_rho_E" + str(0)
    pyplot.plot(y, rho_E_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_E.eps', bbox_inches='tight')
    


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hn", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Process each file in turn
    hdf5files = args

    # FIXME: Remove when plotting multiple files fully supported
    if len(hdf5files)>1:
        print >>sys.stderr, "Only supports 1 file for now.  Sorry."
        return 2
        
    for hdf5file in hdf5files:
        plot(hdf5file)

if __name__ == "__main__":
    sys.exit(main())
