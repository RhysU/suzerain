#!/usr/bin/env python
"""Usage: plot_stats.py HDF5FILE
Plot wall-normal profiles averaged over (X,Z) directions from HDF5FILE.
Options: 
  -f  --file_ext  Output file extension. Default is 'eps'.
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
from scipy.interpolate import interp1d


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

    # Grab reaction rate coefficients
    om_s_coeff = f['bar_om_s'].value
    om_s_coeff = np.array(om_s_coeff).transpose().reshape(Ny,Ns)

    # Grab p coefficients
    p_coeff = f['bar_p'].value
    p_coeff = np.array(p_coeff).reshape(Ny,1)

    # Grab a coefficients
    a_coeff = f['bar_a'].value
    a_coeff = np.array(a_coeff).reshape(Ny,1)

    # Grid delta y, colocation points
    dy = np.diff(y)
    dy = np.append(dy,dy[Ny-2])
    dy = np.array(dy).reshape(Ny,1)

    # Grab mu coefficients
    mu_coeff = f['bar_mu'].value
    mu_coeff = np.array(mu_coeff).reshape(Ny,1)

    # Grab nu coefficients
    nu_coeff = f['bar_nu'].value
    nu_coeff = np.array(nu_coeff).reshape(Ny,1)

    # Grab breakpoints
    yb = f['breakpoints_y'].value

    # Grid delta y, breakpoints
    dyb = np.diff(yb)
    dyb = np.append(dyb,dyb[Ny-6])
    dyb = np.array(dyb).reshape(Ny-4,1)

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
    om_s_col    = D0*om_s_coeff
    mu_col      = D0*mu_coeff
    nu_col      = D0*nu_coeff
    p_col       = D0*p_coeff
    a_col       = D0*a_coeff


    # Computed quantities
    # - Favre averages
    fav_u       = np.array(rho_u_col[:,0]/rho_col[:,0]).reshape(Ny,1)
    fav_v       = np.array(rho_u_col[:,1]/rho_col[:,0]).reshape(Ny,1)
    fav_w       = np.array(rho_u_col[:,2]/rho_col[:,0]).reshape(Ny,1)

    fav_H       = np.array((rho_E_col[:,0] + p_col[:,0])/rho_col[:,0]).reshape(Ny,1)

    # - Bar rho_upp
    # rho_upp     = np.array(np.ravel(rho_col) * np.ravel(fav_u)).reshape(Ny,1)
    rho_upp     = rho_u_col[:,0] - np.array(np.ravel(rho_col) * np.ravel(fav_u)).reshape(Ny,1)
    rho_vpp     = rho_u_col[:,1] - np.array(np.ravel(rho_col) * np.ravel(fav_v)).reshape(Ny,1)
    rho_wpp     = rho_u_col[:,2] - np.array(np.ravel(rho_col) * np.ravel(fav_w)).reshape(Ny,1)

    # - Reynolds stresses
    R_u_u_col   = rho_u_u_col[:,0] - np.multiply(np.multiply(rho_u_col[:,0], rho_u_col[:,0]).reshape(Ny,1), 1/rho_col[:,0])
    R_u_u_col  += 2.0 * np.array(np.ravel(rho_upp) * np.ravel(fav_u)).reshape(Ny,1)

    R_u_v_col   = rho_u_u_col[:,1] - np.multiply(np.multiply(rho_u_col[:,0], rho_u_col[:,1]).reshape(Ny,1), 1/rho_col[:,0])
    R_u_v_col  += np.array(np.ravel(rho_upp) * np.ravel(fav_v)).reshape(Ny,1) + np.array(np.ravel(rho_vpp) * np.ravel(fav_u)).reshape(Ny,1)

    R_v_v_col   = rho_u_u_col[:,3] - np.multiply(np.multiply(rho_u_col[:,1], rho_u_col[:,1]).reshape(Ny,1), 1/rho_col[:,0])
    R_v_v_col  += 2.0 * np.array(np.ravel(rho_vpp) * np.ravel(fav_v)).reshape(Ny,1)

    R_w_w_col   = rho_u_u_col[:,5] - np.multiply(np.multiply(rho_u_col[:,2], rho_u_col[:,2]).reshape(Ny,1), 1/rho_col[:,0])
    R_w_w_col  += 2.0 * np.array(np.ravel(rho_wpp) * np.ravel(fav_w)).reshape(Ny,1)

    # - Viscous ts
    nub = np.interp(yb,y,np.ravel(nu_col))
    nub = np.array(nub).reshape(Ny-4,1)
    inv_nub = 1/nub
    dssqr_over_2nu = np.multiply(dyb[:,0], dyb[:,0]).reshape(Ny-4,1)
    dssqr_over_2nu = np.multiply(dssqr_over_2nu[:,0], 1/nub[:,0]).reshape(Ny-4,1) * 0.5

    # Plots
    figid  = 0
    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_" + str(ifile)
    pyplot.plot(y, rho_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho.' + fileext, bbox_inches='tight')

    figid += 1   
    pyplot.figure(figid)
    key = "bar_rho_u" + str(ifile)
    pyplot.plot(y, rho_u_col[:,0], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_u.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_v" + str(ifile)
    pyplot.semilogx(y, rho_u_col[:,1], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_v.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_w" + str(ifile)
    pyplot.plot(y, rho_u_col[:,2], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_w.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_E" + str(ifile)
    pyplot.plot(y, rho_E_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_E.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_T" + str(ifile)
    pyplot.semilogx(y, T_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_T.' + fileext, bbox_inches='tight')
  
    figid += 1
    pyplot.figure(figid)
    key = "bar_p" + str(ifile)
    pyplot.semilogx(y, p_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_p.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_a" + str(ifile)
    pyplot.semilogx(y, a_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_a.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "R_u_u" + str(ifile)
    pyplot.semilogx(y, R_u_u_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('R_u_u.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "R_u_v" + str(ifile)
    pyplot.semilogx(y, R_u_v_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('R_u_v.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "R_v_v" + str(ifile)
    pyplot.semilogx(y, R_v_v_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('R_v_v.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "R_w_w" + str(ifile)
    pyplot.semilogx(y, R_w_w_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('R_w_w.' + fileext, bbox_inches='tight')

    for s in xrange(0,Ns):
      figid += 1
      pyplot.figure(figid)
      key = "rho_" + sname[s] + "_" + str(ifile)
      pyplot.semilogx(y, rho_s_col[:,s], linewidth=3, label=key)
      pyplot.legend(loc=0)
      rho_file = "rho_" + sname[s] + "." + fileext
      pyplot.savefig(rho_file, bbox_inches='tight')

    for s in xrange(0,Ns):
      figid += 1
      pyplot.figure(figid)
      key = "om_" + sname[s] + "_" + str(ifile)
      pyplot.semilogx(y, om_s_col[:,s], linewidth=3, label=key)
      pyplot.legend(loc=0)
      rho_file = "om_" + sname[s] + "." + fileext
      pyplot.savefig(rho_file, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "dy_collocation" + str(ifile)
    pyplot.loglog(y, dy, 'o-', linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('dy.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "dy_break" + str(ifile)
    pyplot.loglog(yb, dyb, 'o-', linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('dyb.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_mu" + str(ifile)
    pyplot.semilogx(y, mu_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_mu.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_nu" + str(ifile)
    pyplot.semilogx(y, nu_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_nu.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "dssqr_over_2nu" + str(ifile)
    pyplot.loglog(yb, dssqr_over_2nu[:,0], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('dssqr_over_2nu.' + fileext, bbox_inches='tight')

    figid += 1   
    pyplot.figure(figid)
    key = "fav_u" + str(ifile)
    pyplot.plot(y, fav_u[:,0], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('fav_u.' + fileext, bbox_inches='tight')

    figid += 1   
    pyplot.figure(figid)
    key = "fav_H" + str(ifile)
    pyplot.plot(y, fav_H[:,0], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('fav_H.' + fileext, bbox_inches='tight')


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # File extension (eps is default)
    fileext="eps"

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hf:n", ["help", "file_ext="])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
        for o, a in opts:
            # TODO: add sanity check for extension input
	    if o in ("-f", "--file_ext"):
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
