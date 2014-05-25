#!/usr/bin/env python
"""Usage: sgr2suz.py --sgr_file=<sgr_file.dat> --suzerain_file=<suzerain_file.h5>
Interpolate solution from a slow growth RANS simulation to a suzerain 'reacting' file.
Input data must be in text format with the following column layout:
{y, rho_diluter, rho_s, rho_u, rho_v, rho_E}, with rho_diluter the diluter density,
rho_s the species densities (s=1, Ns).
Options:
  -h  --help           This help message.
      --sgr_file=      Text file from sgr.
      --suzerain_file= HDF5 file from suzerain with proper metadata.
"""

import sys
import getopt
import h5py
import numpy as np
import gb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def getsgr2suz(datfile, hdf5file):

    # Load sgr file
    data = np.loadtxt(datfile, delimiter=' ')
    print "Data from sgr loaded"

    # Grab sgr data
    # Number of field variables from sgr data
    nvar   = data.shape[1] - 1

    # Number of grid points in y from sgr
    ny_sgr = data.shape[0]

    # Number of species from sgr
    ns     = nvar - 3  # rho_u, rho_v, rhoE

    # Declare convenience indices for rhou, rhov, rhoE
    irhou = nvar - 2
    irhov = nvar - 1
    irhoE = nvar

    # Declare variables for sgr fields
    y_sgr = np.array(data[:,0    ]).transpose().reshape(ny_sgr,1)
    rho   = np.zeros((ny_sgr,1))
    rho_s = np.zeros((ny_sgr,ns))
    for i in xrange(1,ns+1):
        rho_s[:,i-1] = data[:,i]
        rho         += np.array(data[:,i]).transpose().reshape(ny_sgr,1)

    rho_u = np.array(data[:,irhou]).transpose().reshape(ny_sgr,1)
    rho_v = np.array(data[:,irhov]).transpose().reshape(ny_sgr,1)
    rho_E = np.array(data[:,irhoE]).transpose().reshape(ny_sgr,1)

    # Load suzerain file
    f    = h5py.File(hdf5file,'r+')
    print "Data from suzerain loaded"

    # Grab suzerain mesh info
    Nyf  = f['Ny'].value
    xf   = f['collocation_points_x'].value
    yf   = f['collocation_points_y'].value
    zf   = f['collocation_points_z'].value

    # Grab number of species
    Ns=f['antioch_constitutive_data'].attrs['Ns'][0]

    # Grab relevant metadata
    lower_v =f['lower_v' ].value
    lower_cs=f['lower_cs'].value
    upper_cs=f['upper_cs'].value

    # Grab species names
    sname= np.chararray(Ns, itemsize=5)
    for s in xrange(0,Ns):
      sname[s] = f['antioch_constitutive_data'].attrs['Species_'+str(s)]
      sname[s] = sname[s].strip(' ')

    # Read suzerain data
    frho    = f['rho'].value
    frhou   = f['rho_u'].value
    frhov   = f['rho_v'].value
    frhow   = f['rho_w'].value
    frhoE   = f['rho_E'].value

    # Interpolate data
    frho [:,0,0] = np.interp(yf[:],y_sgr[:,0],rho   [:,0])
    frhou[:,0,0] = np.interp(yf[:],y_sgr[:,0],rho_u [:,0])
    frhov[:,0,0] = np.interp(yf[:],y_sgr[:,0],rho_v [:,0])
    frhoE[:,0,0] = np.interp(yf[:],y_sgr[:,0],rho_E [:,0])

    # Write fields to suzerain file
    f['rho'  ][()]   = frho
    f['rho_u'][()]   = frhou
    f['rho_v'][()]   = frhov
    f['rho_w'][()]   = 0
    f['rho_E'][()]   = frhoE

    # Flow metadata
    lower_v[0]       = frhov[0,0,0]/ frho[0,0,0]
    f['lower_v'][()] = lower_v

    # Species
    frhos   = f['rho'].value  ## Placeholder for species field

    # NOTE: s=0 is the diluter, skip writting it, but compute
    # lower and upper concentrations
    frhos [:,0,0]   = np.interp(yf[:],y_sgr[:,0],rho_s [:,0])
    lower_cs[0]      = frhos[0,0,0]/ frho[0,0,0]
    upper_cs[0]      = frhos[Nyf-1,0,0]/ frho[Nyf-1,0,0]
    for s in xrange(1,Ns):
      rhos_key = 'rho_' + sname[s]
      frhos [:,0,0]   = np.interp(yf[:],y_sgr[:,0],rho_s [:,s])
      f[rhos_key.strip(' ')][()] = frhos
      lower_cs[s]      = frhos[0,0,0]/ frho[0,0,0]
      upper_cs[s]      = frhos[Nyf-1,0,0]/ frho[Nyf-1,0,0]

    # Species metadata
    f['lower_cs'][()] = lower_cs
    f['upper_cs'][()] = upper_cs

    # Close suz file
    f.close()
    print "Data written to suzerain file"


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hn",
                         [ "help"
                           , "sgr_file="
                           , "suzerain_file="
                         ])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            if o in ("--sgr_file"):
                sgr_file = a
                print "sgr file is      : ", sgr_file
            if o in ("--suzerain_file"):
                suzerain_file = a
                print "suzerain file is : ", suzerain_file
        if len(args) > 0:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage as err:
        print >>sys.stderr, err.msg
        return 2

    # Interpolate data form sgr to suz file
    getsgr2suz(sgr_file, suzerain_file)

if __name__ == "__main__":
    sys.exit(main())
