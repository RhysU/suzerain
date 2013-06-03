#!/usr/bin/env python
"""Usage: plot_stats.py HDF5FILE
Plot wall-normal profiles averaged over (X,Z) directions from HDF5FILE.
Options: 
  -f  --file-ext  Output file extension. Default is 'eps'.
  -h  --help      This help message.
"""

# TODO: Add ability to plot things other than bar_rho.

import sys
import getopt
import h5py
import numpy as np
import gb 
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def plot(hdf5file, fileext, ifile):
    print "Plotting", hdf5file

    # Load a stats file
    f = h5py.File(hdf5file,'r')
        
    # Grab number of collocation points and B-spline order
    Ny=f['Ny'].value
    k=f['k'].value
    
    # Grab collocation points
    y = f['collocation_points_y'].value

    # Grab number of species
    Ns=f['antioch_constitutive_data'].attrs['Ns'][0]

    # Grab species names
    sname= np.chararray(Ns, itemsize=5)
    for s in xrange(0,Ns):
      sname[s]=f['antioch_constitutive_data'].attrs['Species_'+str(s)]

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

    # Grab T coefficients
    T_coeff = f['bar_T'].value
    T_coeff = np.array(T_coeff).reshape(Ny,1)

    # Grab rho_u_u coefficients
    rho_u_u_coeff = f['bar_rho_u_u'].value
    rho_u_u_coeff = np.array(rho_u_u_coeff).transpose().reshape(Ny,6)

    # Grab rho_s coefficients
    rho_s_coeff = f['bar_rho_s'].value
    rho_s_coeff = np.array(rho_s_coeff).transpose().reshape(Ny,Ns)

    # Done getting data
    f.close()
    
    D0 = D0T.transpose()

    # Coefficients -> Collocation points
    rho_col     = D0*rho_coeff
    rho_u_col   = D0*rho_u_coeff
    rho_E_col   = D0*rho_E_coeff
    T_col       = D0*T_coeff
    rho_u_u_col = D0*rho_u_u_coeff
    rho_s_col   = D0*rho_s_coeff

    # Computed quantities - Reynolds stresses
    R_u_u_col   = rho_u_u_col[:,0] - np.multiply(np.multiply(rho_u_col[:,0], rho_u_col[:,0]), 1/rho_col[:,0])
    R_v_v_col   = rho_u_u_col[:,3] - np.multiply(np.multiply(rho_u_col[:,1], rho_u_col[:,1]), 1/rho_col[:,0])
    R_w_w_col   = rho_u_u_col[:,5] - np.multiply(np.multiply(rho_u_col[:,2], rho_u_col[:,2]), 1/rho_col[:,0])

    # Plots
    figid  = 0
    figid += 1
    pyplot.figure(figid)
    rho_key = "bar_rho_" + str(ifile)
    pyplot.plot(y, rho_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho.' + fileext, bbox_inches='tight')

    figid += 1   
    pyplot.figure(figid)
    rho_key = "bar_rho_u" + str(ifile)
    pyplot.plot(y, rho_u_col[:,0], linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_u.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    rho_key = "bar_rho_v" + str(ifile)
    pyplot.semilogx(y, rho_u_col[:,1], linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_v.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    rho_key = "bar_rho_w" + str(ifile)
    pyplot.plot(y, rho_u_col[:,2], linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_w.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    rho_key = "bar_rho_E" + str(ifile)
    pyplot.plot(y, rho_E_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_rho_E.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    rho_key = "bar_T" + str(ifile)
    pyplot.plot(y, T_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('bar_T.' + fileext, bbox_inches='tight')
   
    figid += 1
    pyplot.figure(figid)
    rho_key = "R_u_u" + str(ifile)
    pyplot.semilogx(y, R_u_u_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('R_u_u.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    rho_key = "R_v_v" + str(ifile)
    pyplot.semilogx(y, R_v_v_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('R_v_v.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    rho_key = "R_w_w" + str(ifile)
    pyplot.semilogx(y, R_w_w_col, linewidth=3, label=rho_key)
    pyplot.legend(loc=1)
    pyplot.savefig('R_w_w.' + fileext, bbox_inches='tight')

    for s in xrange(0,Ns):
      figid += 1
      pyplot.figure(figid)
      rho_key = "rho_" + sname[s] + "_" + str(ifile)
      pyplot.semilogx(y, rho_s_col[:,s], linewidth=3, label=rho_key)
      pyplot.legend(loc=1)
      rho_file = "rho_" + sname[s] + "." + fileext
      pyplot.savefig(rho_file, bbox_inches='tight')


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # File extension (eps is default)
    fileext="eps"

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hf:n", ["help", "file-ext"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
        for o, a in opts:
            # TODO: add sanity check for input extension
            if o in ("-f", "--file-ext="):
                fileext=a
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Process each file in turn
    hdf5files = args

    # Plot multiple files   
    ifile = 0    
    for hdf5file in hdf5files:
        plot(hdf5file, fileext, ifile)
        ifile += 1

if __name__ == "__main__":
    sys.exit(main())
