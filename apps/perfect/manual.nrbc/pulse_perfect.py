#!/usr/bin/env python
######################################################
# Python script to replace an existing field with
# that of a pulse traveling in the +y direction
# to test boundary conditions. Pulse field taken from
# "Accurate Boundary Conditions for Multicomponent Reacting Flows"
# M. Baum, T. Poinsot, D, Thevenin, J.Comp.Phys. 116, 247-261, 1994
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
f   = h5py.File(sys.argv[1],'r+')

y   = f['collocation_points_y'].value
# print y

# Get these values from file
T_wall = 1
Ly     = f['Ly'].value
Ny     = f['Ny'].value

# Assume these single species parameters are set
gamma = 1.4
R     = 8314.472 / 28.96000
mu    = 1.7e-5

# Mean field parameters
v0    = 0
rho0  = 0.001
p0    = rho0 * R * T_wall
c0    = numpy.sqrt(gamma*R*T_wall)

# Set these values
Ac    = 0.1
Bc    = 10

# Nondimensional
u0    = Ac
Re    = rho0*Ly*u0/mu
Ma    = u0/c0
#print Ma

# Placeholder for auxiliary computed variables
vp    = f['rho_u'].value
pp    = f['rho_u'].value
rhop  = f['rho_u'].value
rhoep = f['rho_u'].value
Tp    = f['rho_u'].value

# v-velocity profile
vp  [:,0,0]  = v0 + Ac * numpy.exp(-numpy.power(Bc/Ly*(y[:]-Ly/2),2))
vwall        = vp[0,0,0]
vp  [:,0,0] -= vwall

# compute pressure, density and energy
pp   [:,0,0] = p0 + rho0 * c0 * (vp[:,0,0]-v0)
rhop [:,0,0] = rho0 + rho0 * (vp[:,0,0]-v0) / c0
rhoep[:,0,0] = pp[:,0,0]/(gamma-1) + numpy.power(vp[:,0,0],2)/2*rhop[:,0,0]

#print rho_u.shape

# Recompute field variables
f['rho'  ][()]  = rhop/rho0
f['rho_u'][()]  = 0
f['rho_v'][()]  = rhop*vp/rho0/u0
f['rho_w'][()]  = 0
f['rho_E'][()]  = rhoep/rho0/c0/c0

# Rescale/reset statistics and other quantities
f['t'         ][()] *= 0
f['bulk_rho'  ][()]  = 1.000884082633707
f['bulk_rho_u'][()] *= 0
f['Re'][()]          = Re
f['Ma'][()]          = Ma

# Close everything
f.close()
