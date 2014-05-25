#!/usr/bin/env python
"""Usage: plot_stats.py HDF5FILE
Plot wall-normal profiles averaged over (X,Z) directions from HDF5FILE.
Options:
  -f  --file_ext  Output file extension. Default is 'eps'.
  -h  --help      This help message.
      --plot_all  Generate secondary 'debug' type of plots.
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


def plot(hdf5file, fileext, ifile, plot_all):
    print "Plotting", hdf5file

    # Load a stats file
    f = h5py.File(hdf5file,'r')

    # Grab number of collocation points and B-spline order
    Ny=f['Ny'].value[0]
    k=f['k'].value[0]

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

    # Get "mass" matrix and convert to dense format
    D1T_gb = f['Dy1T'].value
    D1T = gb.gb2ge(D1T_gb, Ny, k-2)

    # Get "mass" matrix and convert to dense format
    D2T_gb = f['Dy2T'].value
    D2T = gb.gb2ge(D2T_gb, Ny, k-2)

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

    # Grab T coefficients
    T_T_coeff = f['bar_T_T'].value
    T_T_coeff = np.array(T_T_coeff).reshape(Ny,1)

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

    # load baseflow coefficients
    base_rho   = None
    base_rho_u = None
    base_rho_v = None
    base_rho_w = None
    base_rho_E = None
    base_p     = None
    if "largo_baseflow" in f:
      if f['largo_baseflow'].attrs['coefficient_base'] == 'polynomial':
        baseflow_coeff = f['largo_baseflow'].value
        # print 'baseflow coefficients loaded'
        npoly = baseflow_coeff.shape[1]
        base_rho   = np.zeros((Ny,1))
        base_rho_u = np.zeros((Ny,1))
        base_rho_v = np.zeros((Ny,1))
        base_rho_w = np.zeros((Ny,1))
        base_rho_E = np.zeros((Ny,1))
        base_p     = np.zeros((Ny,1))
        for i in xrange(0,npoly):
          for j in xrange(0, Ny):
            y_power_i = np.power(y[j,],i)
            base_rho  [j,0] += y_power_i  * baseflow_coeff[0,i]
            base_rho_u[j,0] += y_power_i  * baseflow_coeff[1,i]
            base_rho_v[j,0] += y_power_i  * baseflow_coeff[2,i]
            base_rho_w[j,0] += y_power_i  * baseflow_coeff[3,i]
            base_rho_E[j,0] += y_power_i  * baseflow_coeff[4,i]
            base_p    [j,0] += y_power_i  * baseflow_coeff[5,i]
      #else:
        # skip loading
        # print 'baseflow coefficients not polynomial'

    # Done getting data
    f.close()

    D0 = D0T.transpose()
    D1 = D1T.transpose()
    D2 = D2T.transpose()

    # Coefficients -> Collocation points
    rho_col     = D0*rho_coeff
    rho_u_col   = D0*rho_u_coeff
    rho_E_col   = D0*rho_E_coeff
    T_col       = D0*T_coeff
    T_T_col     = D0*T_T_coeff
    rho_u_u_col = D0*rho_u_u_coeff
    rho_s_col   = D0*rho_s_coeff
    om_s_col    = D0*om_s_coeff
    mu_col      = D0*mu_coeff
    nu_col      = D0*nu_coeff
    p_col       = D0*p_coeff
    a_col       = D0*a_coeff

    rho_E_col_y  = D1*rho_E_coeff
    rho_E_col_yy = D2*rho_E_coeff
    p_col_y      = D1*p_coeff
    p_col_yy     = D2*p_coeff
    rho_col_y    = D1*rho_coeff
    rho_col_yy   = D2*rho_coeff


    # Computed quantities
    # - Favre averages
    fav_u       = np.array(rho_u_col[:,0]/rho_col[:,0]).reshape(Ny,1)
    fav_v       = np.array(rho_u_col[:,1]/rho_col[:,0]).reshape(Ny,1)
    fav_w       = np.array(rho_u_col[:,2]/rho_col[:,0]).reshape(Ny,1)

    fav_H       = np.array((rho_E_col[:,0] + p_col[:,0])/rho_col[:,0]).reshape(Ny,1)

    if (plot_all):
        # d(d(\fav{H})/dy)/dy
        rho_col_2   =    np.multiply(rho_col  [:,0], rho_col  [:,0]).reshape(Ny,1)
        rho_col_3   =    np.multiply(rho_col_2[:,0], rho_col  [:,0]).reshape(Ny,1)
        rho_col_y2  =    np.multiply(rho_col_y[:,0], rho_col_y[:,0]).reshape(Ny,1)
        fav_H_yy    =    np.array   ((rho_E_col_yy[:,0] + p_col_yy[:,0])/rho_col  [:,0]                 ).reshape(Ny,1)
        fav_H_yy   -= 2.*np.multiply((rho_E_col_y [:,0] + p_col_y [:,0])/rho_col_2[:,0], rho_col_y [:,0]).reshape(Ny,1)
        fav_H_yy   -=    np.multiply((rho_E_col   [:,0] + p_col   [:,0])/rho_col_2[:,0], rho_col_yy[:,0]).reshape(Ny,1)
        fav_H_yy   += 2.*np.multiply((rho_E_col   [:,0] + p_col   [:,0])/rho_col_3[:,0], rho_col_y2[:,0]).reshape(Ny,1)


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

    # - Temperature rms
    Tp_Tp = T_T_col - np.multiply(T_col,T_col).reshape(Ny,1)

    # Plots
    figid  = 0
    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_" + str(ifile)
    if (ifile == 0 and base_rho is not None):
      pyplot.plot(y, base_rho[:,0], linewidth=1)
    pyplot.plot(y, rho_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_u" + str(ifile)
    if (ifile == 0 and base_rho_u is not None):
      pyplot.plot(y, base_rho_u[:,0], linewidth=1)
    pyplot.plot(y, rho_u_col[:,0], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_u.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_v" + str(ifile)
    if (ifile == 0 and base_rho_v is not None):
      pyplot.plot(y, base_rho_v[:,0], linewidth=1)
    pyplot.semilogx(y, rho_u_col[:,1], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_v.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_w" + str(ifile)
    if (ifile == 0 and base_rho_w is not None):
      pyplot.plot(y, base_rho_w[:,0], linewidth=1)
    pyplot.plot(y, rho_u_col[:,2], linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('bar_rho_w.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_rho_E" + str(ifile)
    if (ifile == 0 and base_rho_E is not None):
      pyplot.plot(y, base_rho_E[:,0], linewidth=1)
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
    key = "sqrt_Tp_Tp" + str(ifile)
    pyplot.semilogx(y, np.sqrt(Tp_Tp)/T_col, linewidth=3, label=key)
    pyplot.legend(loc=0)
    pyplot.savefig('sqrt_Tp_Tp.' + fileext, bbox_inches='tight')

    figid += 1
    pyplot.figure(figid)
    key = "bar_p" + str(ifile)
    if (ifile == 0 and base_p is not None):
      pyplot.plot(y, base_p[:,0], linewidth=1)
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

    if (plot_all):
        figid += 1
        pyplot.figure(figid)
        key = "fav_H_yy" + str(ifile)
        pyplot.plot(y, fav_H_yy[:,0], linewidth=0.1, label=key)
        pyplot.axhline(linewidth=0.1, color='r')
        pyplot.legend(loc=0)
        pyplot.savefig('fav_H_yy.' + fileext, bbox_inches='tight')


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # File extension (eps is default)
    fileext="eps"

    # Plot all stuff
    # ... include "debug" type plots one may want to see/declare
    plot_all = False

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hf:n", ["help", "file_ext=", "plot_all"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            if o in ("-f", "--file_ext"):
                fileext=a
            if o in ("-f", "--plot_all"):
                plot_all = True
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage as err:
        print >>sys.stderr, err.msg
        return 2

    # Process each file in turn
    hdf5files = args

    # Plot multiple files
    ifile = 0
    for hdf5file in hdf5files:
        plot(hdf5file, fileext, ifile, plot_all)
        ifile += 1

if __name__ == "__main__":
    sys.exit(main())
