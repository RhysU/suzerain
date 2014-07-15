#!/usr/bin/env python
# Given a job directory compute quick resolution information
# from the last qoi.dat and restart00000.h5 within directory
from __future__ import print_function, division
import h5py
import numpy as np
import sys
import os.path

assert len(sys.argv) == 2
assert os.path.isdir(sys.argv[1])
datf = sys.argv[1]+"/qoi.dat"
hdf5 = sys.argv[1]+"/restart00000.h5"

delta_nu = None
u_tau    = None
delta99  = None
print("qoi.dat:          %s" % datf)
for line in open(datf, 'r'):
    if    line.find('.visc ') != -1:
        s = line.split()
        delta_nu = s[6]
        u_tau    = s[8]
    elif line.find('.thick ') != -1:
        s = line.split()
        delta99 = s[10]

print("delta_nu:         %s" % delta_nu)
delta_nu = float(delta_nu)

print("u_tau:            %s" % u_tau)
u_tau = float(u_tau)

if delta99:
    print("delta99:          %s" % delta99)
    delta99 = float(delta99)
    print("u_tau/delta99:    %s" % (u_tau/delta99))

print("restart00000.h5:  %s" % hdf5)
a = h5py.File(hdf5, 'r')

Nx = a['Nx'][0]
Ny = a['Ny'][0]
Nz = a['Nz'][0]
Lx = a['Lx'][0]
Ly = a['Ly'][0]
Lz = a['Lz'][0]
ht = a['htdelta'][0]
y  = a['collocation_points_y'][()]

print("Nx:               %s" % Nx)
print("Ny:               %s" % Ny)
print("htdelta:          %s" % ht)
print("Nz:               %s" % Nz)

print("\\Delta{}x^{+}:    %s" % (Lx / Nx / delta_nu))
print("\\Delta{}z^{+}:    %s" % (Lz / Nz / delta_nu))

y /= float(delta_nu)
t = np.abs(y - 10).argmin()
print("y_1^{+}:          %s" % (y[ 1] - y[0]))
print("y_{10}^{+}:       %s" % (y[10] - y[0]))
print("y inside 10^{+}:  %s" % t)
print("y_{%d}^{+}:       %s" % (t-1, y[t-1] - y[0]))
print("y_{%d}^{+}:       %s" % (t, y[t] - y[0]))
print("y_{%d}^{+}:       %s" % (t+1, y[t+1] - y[0]))
