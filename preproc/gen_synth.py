import numpy
import h5py
import math
import numpy as np
from scipy.interpolate import griddata

#
# first, generate a synthetic compDNS file
#
g   = h5py.File('synth-comp.h5','w')

st  = ['field_var_RHO','field_var_rhoU','field_var_rhoV','field_var_rhoW','field_var_rhoE']

nx = 10
ny = 20
nz = 30

Lx = 2
Lz = 4

field = np.zeros((nx,ny,nz))

for s in st:
    # create {10,20,30} field of zeros (for a moment)
    g[s] = field

g['mesh_xm'] = np.linspace(0., Lx, nx)
g['mesh_ym'] = np.linspace(-1., 1., ny)
g['mesh_zm'] = np.linspace(0., Lz, nz)

g.close()
    
#
# now, fake suz file
#
f = h5py.File('synth-suz.h5','w')

# populate file with values of problem
f['Lx'] = Lx
f['Nx'] = nx
f['Ly'] = 2
f['Ny'] = ny
f['Lz'] = Lz
f['Nz'] = nz

stt = ['rho','rho_u','rho_v','rho_w','rho_E']
for s in stt:
    # create field of zeros (for a moment)
    f[s] = field

g['collocation_points_x'] = np.linspace(0., Lx, nx)
g['collocation_points_y'] = np.linspace(-1., 1., ny)
g['collocation_points_z'] = np.linspace(0., Lz, nz)

f.close()
#
# nick
# 5/23/13
#
