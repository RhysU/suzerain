######################################################
# Python script to replace an existing field with 
# that of a pulse traveling in the +y direction
# to test boundary conditions
#
# Inputs:
# - The input file name in f (the file will be rewritten)
# - Parameters under 'Set these values': (by default
#   they are set for air)
######################################################

import numpy
import h5py
import math
import numpy as np
from scipy.interpolate import griddata

#
# interpolate from compDNS to suz grid
# and,
# return the interpolated field
#
def interp(c,s, state, field):
    
    flat    = np.array([])

    print 'flatten ', state
    for i in xrange(len(field)):
        for j in xrange(len(field[0])):
            for k in xrange(len(field[0][0])):
                flat    = np.append(flat,field[i][j][k])                
    #
    # two different possible interpolation methods
    # 'nearest' and 'linear'
    #
    print 'interpolating: ', state
    #print len(c), len(flat)
    grid_z0 = griddata(c, flat, s, method='linear')
    return grid_z0

#
# not sure what this function was used for
#
def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

######################################################
# Open all of the files
######################################################
#
# suzerain file
#
#f   = h5py.File('temporal_M0p3_physical.h5','r+')
f   = h5py.File('synth-suz.h5','r+')
#
# compDNS file
#
g   = h5py.File('synth-comp.h5','r+')
#g   = h5py.File('data_6.8800011442E-01.h5','r+')

## Get these values from file (account for 3/2 dealiasing)
Lx     = f['Lx'].value
Nx     = f['Nx'].value * 3 / 2.
Ly     = f['Ly'].value
Ny     = f['Ny'].value
Lz     = f['Lz'].value
Nz     = f['Nz'].value * 3 /2.

## suz variables
frho  = f['rho'].value
frhou = f['rho_u'].value
frhov = f['rho_v'].value
frhow = f['rho_w'].value
frhoE = f['rho_E'].value

## comp variables
grho  = g['field_var_RHO'].value
grhou = g['field_var_rhoU'].value
grhov = g['field_var_rhoV'].value
grhow = g['field_var_rhoW'].value
grhoE = g['field_var_rhoE'].value

#
# Following method of: 
# http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
#

# this is current grid (compDNS)
xg   = g['mesh_xm'].value
yg   = g['mesh_ym'].value
zg   = g['mesh_zm'].value
print "Interpolating  :            ", len(zg),len(yg),len(xg)

# this is the final grid (suz)
grid_x = f['collocation_points_x'].value
grid_z = f['collocation_points_z'].value
grid_y = f['collocation_points_y'].value
print "To grid size of:            ", len(grid_z),len(grid_y),len(grid_x)
print len(grho),len(grho[0]),len(grho[0][0])

# CompDNS field stored as: (z,y,x)
# and
# Suz field stored as    : (y,z,x)
#
cpoints = np.array([np.array([xg[0], yg[0], zg[0]])])
spoints = np.array([np.array([grid_x[0], grid_y[0], grid_z[0]])])
flat    = np.array([])

#
# compDNS grid -- flatten 
#
print 'flatten compDNS'
for i in xrange(len(xg)):
    print i 
    for j in xrange(len(yg)):
        for k in xrange(len(zg)):
            cpoints = np.vstack((cpoints,np.array([xg[i], yg[j], zg[k]])))

#
# flatten suzerain grid 
# (this is the grid we are interpolating to)
#
print 'flatten Suz'
for i in xrange(len(grid_x)):
    print i
    for j in xrange(len(grid_y)):
        for k in xrange(len(grid_z)):
            spoints = np.vstack((spoints,np.array([grid_x[i], grid_y[j], grid_z[k]])))

# hack because numpy wont let me cat two uneven arrays
cpoints = np.delete(cpoints, 0, 0)
spoints = np.delete(spoints, 0, 0)

# iterate over each state variable we need to interpolate
st  = ['field_var_RHO','field_var_rhoU','field_var_rhoV','field_var_rhoW','field_var_rhoE']
stt = ['rho','rho_u','rho_v','rho_w','rho_E']

# counter
ct = 0

for s in st:
    # grab state var
    field = g[s].value

    # interpolate
    interp_g = interp(cpoints, spoints, s, field)

    # save interpolated grid to suz
    #print stt[ct], len(interp_g)
    #print interp_g
    f[stt[ct]] = interp_g

    # iterate counter
    ct = ct + 1

# old method: 
#
# # Map comp field to suzerain
# # TODO: interpolate results
# # currently, it is just a direct mapping. 
# for k in xrange(0,Nz):
#   for j in xrange(0,Ny):
#     for i in xrange(0,Nx):
#       #  print i,j,k
#       frho [j,k,i] = grho [k,j,i]
#       frhou[j,k,i] = grhou[k,j,i]
#       frhov[j,k,i] = grhov[k,j,i]
#       frhow[j,k,i] = grhow[k,j,i]
#       frhoE[j,k,i] = grhoE[k,j,i]

# TODO: comp stores the values at cell centers,
#       therefore the values at the wall boundaries
#       are not in the fields. 
#       Maybe reset the wall boundary values 
#       (I think this will be done at run time by suz in any case)

# Print streamwise momentum pencil
#print frhou[:,0,0]

# # Assign interpolated fields to suz
# f['rho'  ][()] = frho 
# f['rho_u'][()] = frhou
# f['rho_v'][()] = frhov
# f['rho_w'][()] = frhow
# f['rho_E'][()] = frhoE

# Close everything
f.close()
g.close()
print "DONE!"

#
# written by victor
# heavily modified by nick
# 
# 5/22/13
#
