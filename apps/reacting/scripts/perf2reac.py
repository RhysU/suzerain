######################################################
# Python script to convert a nondimensional channel field
# from the 'perfect' implementation to a dimensional field
# for the 'reacting' implementation
# 
# Inputs:
# - The input file name in f (the file will be rewritten)
# - Parameters under 'Set these values': (by default
#   they are set for air)
######################################################

import numpy
import h5py
import math
 
######################################################
# Open all of the files
f   = h5py.File('physical831181_0_reacting.h5','r+')

# Relevant parameters from perfect
Re     = f['/Re'].value
Ma     = f['/Ma'].value
gamma  = f['/gamma'].value
beta   = f['/beta'].value
Tw     = 1.

######################################################
# NOTE: Set these values
#
# Cp
Cp_d = 1.0045e+3

# Dimensional reference temperature and viscosity (for air)
T0_d     = 273
mu0_d    = 1.716e-5

# Dimensional wall temperature
T_wall_d = 273

#####################################################
# Dimensional viscosity at T_wall
mu0d_w    = mu0_d * numpy.power(T_wall_d/T0_d, beta)

# Relevant input dimensional parameters
Cp         = f.create_dataset('Cp'         , shape=(1,), dtype=('float64'), data=Cp_d)
aux        = Cp.value / gamma
Cv         = f.create_dataset('Cv'         , shape=(1,), dtype=('float64'), data=aux)
T0         = f.create_dataset('T0'         , shape=(1,), dtype=('float64'), data=T0_d)
mu0        = f.create_dataset('mu0'        , shape=(1,), dtype=('float64'), data=mu0_d)
aux        = T_wall_d * Tw
T_wall     = f.create_dataset('T_wall'     , shape=(1,), dtype=('float64'), data=aux)
filter_phi = f.create_dataset('filter_phi' , shape=(1,), dtype=('float64'), data=0)

# Chosen dimensional lenghtscale
l0 = 1

# Computed/dependent dimensional parameters
R      = Cp.value - Cv.value
a0     = numpy.sqrt(gamma*R*T_wall.value)
u0     = Ma * a0
rho0   = Re*mu0d_w/u0/l0
t0     = l0/u0
E0     = a0*a0

# Rescale field variables
f['rho'  ][()] *= rho0
f['rho_u'][()] *= rho0*u0
f['rho_v'][()] *= rho0*u0
f['rho_w'][()] *= rho0*u0
f['rho_E'][()] *= rho0*E0

# Rescale/reset statistics and other quantities
f['t'         ][()] *= 0 #t0
f['bulk_rho'  ][()] *= rho0
f['bulk_rho_u'][()] *= rho0*u0

# Delete unused parameters from the perfect case
del f['/Re']
del f['/Ma']
del f['/gamma']

# Close everything
f.close()
