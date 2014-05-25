#!/usr/bin/env python
"""Usage: plan_suzerain.py HDF5FILE
Get DNS planing information from a slow growth RANS solution.
Options:
  -h  --help                            This help message.
      --dt                              Estimated time-step
      --Lx_over_delta_tgt               Target Lx over boundary layer thickness (default is 10)
      --Ly_over_delta_tgt               Target Ly over boundary layer thickness (default is 2.5)
      --Lz_over_delta_tgt               Target Lz over boundary layer thickness (default is 3)
      --Dx_over_delta_nu_tgt            Target Dx^+ (Dx over viscous length scale, default is 14)
      --Dz_over_delta_nu_tgt            Target Dz^+ (Dz over viscous length scale, default is 7)
      --CPUh_performance                Cluster performance (CPU hours per million points per thousand time steps, default is 170)
      --delta_factor                    Factor to scale the boundary layer thickness (default is 1)
      --y1_plus_tgt                     Target location of first break point in plus units (default is 0.6)
      --Dy_edge_over_delta_ref_tgt      Target mesh refinement at the boundary layer edge (Dy over reference thickness, default is 0.15 of displacement thickness)
                                        NOTE: former option --Dy_edge_over_delta_star_tgt remains valid for backwards compatibility
      --use_theta_reference             Use momentum thickness as reference thickness as outer resolution measure (default is False)
      --samples_per_turnover            Target number of samples per turnover (default is 10)

"""

import sys
import getopt
import h5py
import numpy as np
import gb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def gety(b, L, x):
    y = L * (1 + np.tanh(b*(x - 1))/np.tanh(b))
    return y


def getplan(hdf5file, dt):
    print "Processing", hdf5file

    # Load a stats file
    f = h5py.File(hdf5file,'r')
    # print "File loaded"

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

    # Grab rho_u_u coefficients
    rho_u_u_coeff = f['bar_rho_u_u'].value
    rho_u_u_coeff = np.array(rho_u_u_coeff).transpose().reshape(Ny,6)

    # Grab rho_s coefficients
    rho_s_coeff = f['bar_rho_s'].value
    rho_s_coeff = np.array(rho_s_coeff).transpose().reshape(Ny,Ns)

    # Grab bar_u coefficients
    bar_u_coeff = f['bar_u'].value
    bar_u_coeff = np.array(bar_u_coeff).transpose().reshape(Ny,3)

    # Grab bar_u coefficients
    bar_mu_coeff = f['bar_mu'].value
    bar_mu_coeff = np.array(bar_mu_coeff).transpose().reshape(Ny,1)

    # Grab bar_u coefficients
    bar_a_coeff = f['bar_a'].value
    bar_a_coeff = np.array(bar_a_coeff).transpose().reshape(Ny,1)

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
    rho_u_u_col = D0*rho_u_u_coeff
    rho_s_col   = D0*rho_s_coeff
    bar_u_col   = D0*bar_u_coeff
    bar_du_col  = D1*bar_u_coeff

    # Input parameters
    T_wall     =  T_coeff[0,0]
    T_inf      =  T_coeff[Ny-1,0]
    V_wall     =  bar_u_coeff[0,1]
    U_inf      =  bar_u_coeff[Ny-1,0]
    mu_wall    =  bar_mu_coeff[0,0]
    mu_inf     =  bar_mu_coeff[Ny-1,0]
    a_inf      =  bar_a_coeff[Ny-1,0]
    rho_u_inf  =  rho_u_coeff[Ny-1,0]

    # Raw parameters
    rho_wall   =  rho_coeff[0,0]
    rho_inf    =  rho_coeff[Ny-1,0]
    dudy_wall  =  bar_du_col[0,0]


    # Compute bl thickness
    # NOTE: keeping the value from comp for now,
    # since the simulation is not converged
    jdelta = 0
    for j in xrange(0,Ny):
        if (bar_u_col[j,0] < 0.99*U_inf):
           jdelta += 1
        else:
           break
    delta = y[jdelta]

    # FIXME: Generalize definition of edge values
    rho_edge  = rho_inf
    rhoU_edge = rho_u_inf
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

    # Pick thickness as reference for edge resolution
    if use_theta_reference :
        delta_ref_resolution = theta
    else:
        delta_ref_resolution = delta_star

    # Computed parameters
    Re_delta_star = rho_inf * U_inf * delta_star / mu_inf
    Re_theta      = rho_inf * U_inf * theta      / mu_inf
    H1            = delta_star / theta
    H2            = delta      / theta
    tau_wall      = mu_wall * dudy_wall
    u_tau         = np.sqrt(tau_wall / rho_wall)
    delta_nu      = mu_wall / rho_wall / u_tau
    Re_tau        = rho_wall * u_tau * delta / mu_wall
    turnover_time = delta / u_tau

    # Output parameters
    print 'Boundary layer input parameters'
    print 'T_wall                      = ', T_wall
    print 'T_inf                       = ', T_inf
    print 'U_inf                       = ', U_inf
    print 'V_wall                      = ', V_wall
    print
    print 'Raw parameters'
    print 'Ma_inf                      = ', U_inf/a_inf
    print 'rho_wall                    = ', rho_wall
    print 'rho_inf                     = ', rho_inf
    print 'delta_star                  = ', delta_star
    print 'theta                       = ', theta
    print 'delta                       = ', delta
    if dt != 0:
        print 'dt                          = ', dt

    print
    print 'Computed parameters'
    print 'Re_delta_star               = ', Re_delta_star
    print 'Re_theta                    = ', Re_theta
    print 'H1                          = ', H1
    print 'H2                          = ', H2
    print 'dU/dy|_wall                 = ', dudy_wall
    print 'tau|_wall                   = ', tau_wall
    print 'u_tau                       = ', u_tau
    print 'delta_nu                    = ', delta_nu
    print 'Re_tau                      = ', Re_tau
    print 'turnover_time               = ', turnover_time

    # Compute run parameters
    delta_tgt = delta * delta_factor
    Lx_tgt = Lx_over_delta_tgt * delta_tgt
    Ly_tgt = Ly_over_delta_tgt * delta_tgt
    Lz_tgt = Lz_over_delta_tgt * delta_tgt
    Dx_tgt = Dx_over_delta_nu_tgt * delta_nu
    Dz_tgt = Dz_over_delta_nu_tgt * delta_nu

    # Compute y-mesh parameters
    nymin = 0;
    nymax = 4096;
    ny_mid = (nymin+nymax)/2;
    tol_ny = 0.01;
    err_ny = 2*tol_ny;
    while (abs(err_ny) > tol_ny):
      rmin = 0.1
      rmax = 100
      r_mid = (rmin+rmax)/2
      tol_r = 0.01
      err_r = 2*tol_r
      # r is the stretching factor;
      # solve for the r that satisfies the y1_plus condition
      # for this
      while (abs(err_r) > tol_r and (rmax - rmin) > 1.0E-10):
        dy1 = 1/float(ny_mid)
        y_grid = np.zeros((ny_mid,1)).reshape(ny_mid,1)
        for j in xrange(0,ny_mid):
          xloc        = j * dy1
          y_grid[j,0] = gety(r_mid, Ly_tgt, xloc)

        Dy_vec = np.zeros((ny_mid,1)).reshape(ny_mid,1)
        for j in xrange(0,ny_mid-1):
          Dy_vec[j,0] = y_grid[j+1,0]-y_grid[j,0]

        y1_plus_mid = y_grid[1,0]/delta_nu

        err_r = (y1_plus_mid-y1_plus_tgt)/y1_plus_tgt
        if (err_r > tol_r):
           rmin = r_mid
           r_mid = (rmax+r_mid)/2
        elif (err_r < tol_r):
           rmax = r_mid;
           r_mid = (rmin+r_mid)/2
        #else:
          # result is converged

      # Assign the solution to r
      r = r_mid

      # now evaluate if for this grid the criterion for dymax is satisfied
      for j in xrange(0,ny_mid):
        if (y_grid[j,0] < delta_tgt):
          Dy_edge_over_delta_ref_mid = Dy_vec[j,0] / delta_ref_resolution
        else:
          break

      Dy_max_over_delta_ref_mid = (max(Dy_vec) / delta_star)[0]
      err_ny = (Dy_edge_over_delta_ref_mid-Dy_edge_over_delta_ref_tgt)/Dy_edge_over_delta_ref_mid
      # if the dy at edge is larger than the target,
      # increase the number of points,
      # decrease otherwise
      if (err_ny > tol_ny):
         nymin = ny_mid
         ny_mid = (nymax+ny_mid)/2
      elif (err_ny < tol_ny):
         nymax = ny_mid
         ny_mid = (nymin+ny_mid)/2
      #else:
        # result is converged

    # Assign grid parameters when converged
    Ny_tgt  = ny_mid
    htdelta = r_mid

    # Compute target run parameters
    Nx_tgt = Lx_tgt / Dx_tgt
    Nz_tgt = Lz_tgt / Dz_tgt
    Np     = Nx_tgt * Ny_tgt * Nz_tgt
    turnover_tgt = turnover_time * delta_factor

    Nt_turntime       = 0
    CPUh_turnover     = 0
    Nt_turntime_tgt   = 0
    CPUh_turnover_tgt = 0
    if dt != 0:
        Nt_turntime_tgt   = turnover_tgt / dt
        CPUh_turnover_tgt = CPUh_performance * Np * Nt_turntime_tgt / 1.0e9

    dt_sample         = turnover_tgt / samples_per_turnover
    nt_sample         = dt_sample / dt

    # Output run parameters
    print
    print 'Run target parameters'
    print 'delta_factor                = ', delta_factor
    print 'turnover_tgt                = ', turnover_tgt
    print 'Nt_turnover_tgt             = ', Nt_turntime_tgt
    print 'Lx_over_delta_tgt           = ', Lx_over_delta_tgt
    print 'Ly_over_delta_tgt           = ', Ly_over_delta_tgt
    print 'Lz_over_delta_tgt           = ', Lz_over_delta_tgt
    print 'Dx_over_delta_nu_tgt        = ', Dx_over_delta_nu_tgt
    print 'Dz_over_delta_nu_tgt        = ', Dz_over_delta_nu_tgt
    print 'y1_plus_tgt                 = ', y1_plus_tgt
    if (use_theta_reference):
        print 'Reference outer thickness   = ', 'theta'
    else:
        print 'Reference outer thickness   = ', 'delta_star'
    print 'Dy_edge_over_delta_ref_tgt  = ', Dy_edge_over_delta_ref_tgt
    print
    print 'Run setup and cost parameters'
    print 'htdelta                     = ', htdelta
    print 'Lx_tgt                      = ', Lx_tgt
    print 'Ly_tgt                      = ', Ly_tgt
    print 'Lz_tgt                      = ', Lz_tgt
    print 'Nx_tgt                      = ', Nx_tgt
    print 'Ny_tgt                      = ', Ny_tgt
    print 'Nz_tgt                      = ', Nz_tgt
    print 'y1_plus                     = ', y1_plus_mid
    print 'Dy_edge_over_delta_ref      = ', Dy_edge_over_delta_ref_mid
    print 'Dy_max_over_delta_ref       = ', Dy_max_over_delta_ref_mid
    print 'Np_million                  = ', Np / 1.0e6
    print 'CPUh/Np(10^6)/Nt(10^3)      = ', CPUh_performance
    print 'CPUh/turnover               = ', CPUh_turnover_tgt
    print 'samples_per_turnover        = ', samples_per_turnover
    print 'dt_sample                   = ', dt_sample
    print 'nt_sample                   = ', nt_sample



def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Default optional parameters
    global Lx_over_delta_tgt
    global Ly_over_delta_tgt
    global Lz_over_delta_tgt
    global Dx_over_delta_nu_tgt
    global Dz_over_delta_nu_tgt
    global CPUh_performance
    global delta_factor
    global y1_plus_tgt
    global Dy_edge_over_delta_ref_tgt
    global samples_per_turnover
    global use_theta_reference

    dt                          =    0
    Lx_over_delta_tgt           =   10.0
    Ly_over_delta_tgt           =    2.5
    Lz_over_delta_tgt           =    3.0
    Dx_over_delta_nu_tgt        =   14
    Dz_over_delta_nu_tgt        =    7
    CPUh_performance            =  170
    delta_factor                =    1.0
    y1_plus_tgt                 =    0.6
    Dy_edge_over_delta_ref_tgt  =    0.15
    samples_per_turnover        =   10
    use_theta_reference         =    False

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hn",
                         [ "help"
                         , "dt="
                         , "Lx_over_delta_tgt="
                         , "Ly_over_delta_tgt="
                         , "Lz_over_delta_tgt="
                         , "Dx_over_delta_nu_tgt="
                         , "Dz_over_delta_nu_tgt="
                         , "CPUh_performance="
                         , "delta_factor="
                         , "y1_plus_tgt="
                         , "Dy_edge_over_delta_star_tgt="
                         , "Dy_edge_over_delta_ref_tgt="
                         , "samples_per_turnover="
                         , "use_theta_reference"
                         ])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            if o in ("--dt"):
                dt = float(a)
            if o in ("--Lx_over_delta_tgt"):
                Lx_over_delta_tgt = float(a)
            if o in ("--Ly_over_delta_tgt"):
                Ly_over_delta_tgt = float(a)
            if o in ("--Lz_over_delta_tgt"):
                Lz_over_delta_tgt = float(a)
            if o in ("--Dx_over_delta_nu_tgt"):
                Dx_over_delta_nu_tgt = float(a)
            if o in ("--Dz_over_delta_nu_tgt"):
                Dz_over_delta_nu_tgt = float(a)
            if o in ("--CPUh_performance"):
                CPUh_performance = float(a)
            if o in ("--delta_factor"):
                delta_factor = float(a)
            if o in ("--y1_plus_tgt"):
                y1_plus_tgt = float(a)
            # keep "delta_star" option for backwards-compatibility
            if (o in ("--Dy_edge_over_delta_star_tgt") or
                o in ("--Dy_edge_over_delta_ref_tgt" )   ):
                Dy_edge_over_delta_ref_tgt = float(a)
            if o in ("--samples_per_turnover"):
                samples_per_turnover = float(a)
            if o in ("--use_theta_reference"):
                use_theta_reference = True
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Process each file in turn
    hdf5files = args

    # FIXME: Remove when processing multiple files fully supported
    if len(hdf5files)>1:
        print >>sys.stderr, "Only supports 1 file for now.  Sorry."
        return 2

    # Get averaged parameters from multiple files
    for hdf5file in hdf5files:
        getplan(hdf5file, dt)

if __name__ == "__main__":
    sys.exit(main())
