#!/usr/bin/env python
"""Usage: blparam_stats.py HDF5FILE
Get boundary layer parameters from HDF5FILE, and generate output
file (blparam_vs_time) and plots (BL_*).
Options:
  -f  --file_ext  Output file extension. Default is 'eps'.
  -h  --help      This help message.
"""

import sys
import getopt
import h5py
import numpy as np
import gb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def getblparam(hdf5file):
    print "Processing", hdf5file

    # Load a stats file
    f = h5py.File(hdf5file,'r')

    # Grab time
    t=f['t'].value[0]

    # Grab number of points in x and z
    Nx=f['Nx'].value[0]
    Nz=f['Nz'].value[0]

    # Grab domain lengths
    Lx=f['Lx'].value[0]
    Ly=f['Ly'].value[0]
    Lz=f['Lz'].value[0]

    # Grab number of collocation points and B-spline order
    Ny=f['Ny'].value[0]
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

    # Get "mass" matrix and convert to dense format
    D1T_gb = f['Dy1T'].value
    D1T = gb.gb2ge(D1T_gb, Ny, k-2)

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

    # Grab rho_s coefficients
    rho_s_coeff = f['bar_rho_s'].value
    rho_s_coeff = np.array(rho_s_coeff).transpose().reshape(Ny,Ns)

    # Grab bar_u coefficients
    bar_u_coeff = f['bar_u'].value
    bar_u_coeff = np.array(bar_u_coeff).transpose().reshape(Ny,3)

    # Grab bar_u coefficients
    bar_mu_coeff = f['bar_mu'].value
    bar_mu_coeff = np.array(bar_mu_coeff).transpose().reshape(Ny,1)

    # Grab p coefficients
    p_coeff = f['bar_p'].value
    p_coeff = np.array(p_coeff).reshape(Ny,1)

    # Grab a coefficients
    a_coeff = f['bar_a'].value
    a_coeff = np.array(a_coeff).reshape(Ny,1)

    # Grid delta y, collocation
    dy = np.diff(y)
    dy = np.append(dy,dy[Ny-2])
    dy = np.array(dy).reshape(Ny,1)

    # Grab breakpoints
    yb = f['breakpoints_y'].value

    # Grid delta y, breakpoints
    dyb = np.diff(yb)
    dyb = np.append(dyb,dyb[Ny-6])
    dyb = np.array(dyb).reshape(Ny-4,1)

    # Integration weights
    i_weights = f['integration_weights'].value
    i_weights = np.array(i_weights).reshape(Ny,1)

    # Done getting data
    f.close()

    D0 = D0T.transpose()
    D1 = D1T.transpose()

    # Coefficients -> Collocation points
    rho_col     = D0*rho_coeff
    rho_u_col   = D0*rho_u_coeff
    rho_E_col   = D0*rho_E_coeff
    T_col       = D0*T_coeff
    p_col       = D0*p_coeff
    a_col       = D0*a_coeff
    rho_s_col   = D0*rho_s_coeff
    bar_u_col   = D0*bar_u_coeff
    bar_du_col  = D1*bar_u_coeff

    # Input parameters
    T_wall     =  T_coeff[0,0]
    T_inf      =  T_coeff[Ny-1,0]
    V_wall     =  bar_u_coeff[0,1]
    U_inf      =  bar_u_coeff[Ny-1,0]
    V_inf      =  bar_u_coeff[Ny-1,1]
    W_inf      =  bar_u_coeff[Ny-1,2]

    # Grid parameters
    y1         =  y[1]
    y1b        =  yb[1]

    # Raw parameters
    p_wall     =  p_coeff[0,0]
    p_inf      =  p_coeff[Ny-1,0]
    a_wall     =  a_coeff[0,0]
    a_inf      =  a_coeff[Ny-1,0]
    rho_wall   =  rho_coeff[0,0]
    mu_wall    =  bar_mu_coeff[0,0]
    rho_inf    =  rho_coeff[Ny-1,0]
    mu_inf     =  bar_mu_coeff[Ny-1,0]
    dudy_wall  =  bar_du_col[0,0]

    # Compute delta (BL thickness)
    jdelta = 0
    for j in xrange(0,Ny):
        if (bar_u_col[j,0] < 0.99*U_inf):
           jdelta += 1
        else:
           break

    # compute delta by interpolation
    frac  = (0.99*U_inf - bar_u_col[j,0]) / (bar_u_col[jdelta+1,0] - bar_u_col[j,0])
    delta = (y[jdelta+1] - y[jdelta]) * frac + y[jdelta]

    # TODO: Here the value at the edge is defined as the value at "infinity"
    #       This needs to be generalized for baseflow computations
    rho_edge  = rho_inf
    rhoU_edge = rho_inf * U_inf
    U_edge    = U_inf
    mu_edge   = mu_inf

    # Compute delta_star (displacement thickness)
    delta_star = 0
    for j in xrange(0,Ny):
        delta_star += (1 - rho_u_col[j,0] / rhoU_edge) * i_weights[j,0]

    # Compute theta (momentum thickness)
    theta = 0
    for j in xrange(0,Ny):
        theta += rho_u_col[j,0] / rhoU_edge * (1 - bar_u_col[j,0] / U_edge) * i_weights[j,0]

    # Computed parameters
    Re_delta_star    = rho_edge * U_edge * delta_star / mu_edge
    Re_theta         = rho_edge * U_edge * theta      / mu_edge
    Re_delta         = rho_edge * U_edge * delta      / mu_edge
    H1               = delta_star / theta
    H2               = delta      / theta
    tau_wall         = mu_wall * dudy_wall
    u_tau            = np.sqrt(tau_wall / rho_wall)
    delta_nu         = mu_wall / rho_wall / u_tau
    Re_tau           = rho_wall * u_tau * delta / mu_wall
    turnover_time    = delta / u_tau
    flowthrough_time = Lx / U_edge

    Lx_over_delta    = Lx / delta
    Ly_over_delta    = Ly / delta
    Lz_over_delta    = Lz / delta
    y1_plus          = y1 / delta_nu
    y1b_plus         = y1b / delta_nu
    Dx               = Lx / Nx
    Dz               = Lz / Nz
    Dx_plus          = Dx / delta_nu
    Dz_plus          = Dz / delta_nu

    # Resolution in y, collocation
    Ny_below_5plus  = 0
    for j in xrange(0,Ny):
        if (y[j] < 5*delta_nu):
           Ny_below_5plus += 1
        else:
           break

    Ny_below_10plus = 0
    for j in xrange(0,Ny):
        if (y[j] < 10*delta_nu):
           Ny_below_10plus += 1
        else:
           break

    Ny_below_delta = 0
    for j in xrange(0,Ny):
        if (y[j] < delta):
           Ny_below_delta += 1
        else:
           break

    # Resolution in y, breakpoints
    Nyb_below_5plus  = 0
    for j in xrange(0,Ny):
        if (yb[j] < 5*delta_nu):
           Nyb_below_5plus += 1
        else:
           break

    Nyb_below_10plus = 0
    for j in xrange(0,Ny):
        if (yb[j] < 10*delta_nu):
           Nyb_below_10plus += 1
        else:
           break

    Nyb_below_delta = 0
    for j in xrange(0,Ny):
        if (yb[j] < delta):
           Nyb_below_delta += 1
        else:
           break

    # Put parameters in an array
    prms = np.empty([0,1])
    prms = np.append(prms, [t                 ])
    prms = np.append(prms, [rho_wall          ])
    prms = np.append(prms, [p_wall            ])
    prms = np.append(prms, [a_wall            ])
    prms = np.append(prms, [mu_wall           ])
    prms = np.append(prms, [mu_inf            ])
    prms = np.append(prms, [delta_star        ])
    prms = np.append(prms, [theta             ])
    prms = np.append(prms, [delta             ])
    prms = np.append(prms, [H1                ])
    prms = np.append(prms, [H2                ])
    prms = np.append(prms, [dudy_wall         ])
    prms = np.append(prms, [tau_wall          ])
    prms = np.append(prms, [u_tau             ])
    prms = np.append(prms, [delta_nu          ])
    prms = np.append(prms, [y1b_plus          ])
    prms = np.append(prms, [Re_tau            ])
    prms = np.append(prms, [Re_delta_star     ])
    prms = np.append(prms, [Re_theta          ])
    prms = np.append(prms, [Re_delta          ])
    prms = np.append(prms, [U_edge            ])
    prms = np.append(prms, [V_inf             ])
    prms = np.append(prms, [W_inf             ])
    prms = np.append(prms, [T_inf             ])
    prms = np.append(prms, [rho_inf           ])
    prms = np.append(prms, [p_inf             ])
    prms = np.append(prms, [a_inf             ])
    prms = np.append(prms, [Dx_plus           ])
    prms = np.append(prms, [Dz_plus           ])
    prms = np.append(prms, [Lx_over_delta     ])
    prms = np.append(prms, [Ly_over_delta     ])
    prms = np.append(prms, [Lz_over_delta     ])
    prms = np.append(prms, [turnover_time     ])
    prms = np.append(prms, [flowthrough_time  ])

    return prms


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
        except getopt.error as msg:
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
    except Usage as err:
        print >>sys.stderr, err.msg
        return 2

    # Ouptut filename
    outfile = 'blparam_vs_time.out'

    # Create and save header
    head_table = np.empty([0])
    head_table = np.append(head_table, ["t                       "])
    head_table = np.append(head_table, ["rho_wall                "])
    head_table = np.append(head_table, ["p_wall                  "])
    head_table = np.append(head_table, ["a_wall                  "])
    head_table = np.append(head_table, ["mu_wall                 "])
    head_table = np.append(head_table, ["mu_inf                  "])
    head_table = np.append(head_table, ["delta_star              "])
    head_table = np.append(head_table, ["theta                   "])
    head_table = np.append(head_table, ["delta                   "])
    head_table = np.append(head_table, ["H1                      "])
    head_table = np.append(head_table, ["H2                      "])
    head_table = np.append(head_table, ["dudy_wall               "])
    head_table = np.append(head_table, ["tau_wall                "])
    head_table = np.append(head_table, ["u_tau                   "])
    head_table = np.append(head_table, ["delta_nu                "])
    head_table = np.append(head_table, ["y1b_plus                "])
    head_table = np.append(head_table, ["Re_tau                  "])
    head_table = np.append(head_table, ["Re_delta_star           "])
    head_table = np.append(head_table, ["Re_theta                "])
    head_table = np.append(head_table, ["Re_delta                "])
    head_table = np.append(head_table, ["U_inf                   "])
    head_table = np.append(head_table, ["V_inf                   "])
    head_table = np.append(head_table, ["W_inf                   "])
    head_table = np.append(head_table, ["T_inf                   "])
    head_table = np.append(head_table, ["rho_inf                 "])
    head_table = np.append(head_table, ["p_inf                   "])
    head_table = np.append(head_table, ["a_inf                   "])
    head_table = np.append(head_table, ["Dx_plus                 "])
    head_table = np.append(head_table, ["Dz_plus                 "])
    head_table = np.append(head_table, ["Lx_over_delta           "])
    head_table = np.append(head_table, ["Ly_over_delta           "])
    head_table = np.append(head_table, ["Lz_over_delta           "])
    head_table = np.append(head_table, ["turnover_time           "])
    head_table = np.append(head_table, ["flowthrough_time        "])
    nprms      = head_table.shape[0]
    head_table = np.array(head_table).reshape(1, nprms)
    np.savetxt(outfile, head_table, delimiter=' ', fmt='%s')

    # Process each file in turn
    hdf5files = args
    nfiles    = len(hdf5files)

    # Get averaged parameters from multiple files
    # and append the results to the array
    prms_table = np.empty([0])
    for hdf5file in hdf5files:
        prms = getblparam(hdf5file)
        prms_table = np.append(prms_table, prms, axis=0)

    # Reshape table and sort by time
    prms_table = np.array(prms_table).reshape(nfiles, nprms)
    prms_table = prms_table[prms_table[:,0].argsort()]

    # Save table to file
    f_handle = file(outfile, 'a')
    np.savetxt(f_handle, prms_table, delimiter=' ')
    f_handle.close()

    # Generate plots
    print "Generating plots"
    for figid in xrange(1,nprms):
        pyplot.figure(figid)
        key = head_table[0,figid].strip()
        pyplot.plot(prms_table[:,0], prms_table[:,figid], '-o',linewidth=3, label=key)
        pyplot.legend(loc=0)
        pyplot.savefig('BL_' + key + '.' + fileext, bbox_inches='tight')

if __name__ == "__main__":
    sys.exit(main())
