#!/usr/bin/env python
from __future__ import print_function, division
import h5py
import numpy as np
import sys

assert len(sys.argv) == 3

a        = h5py.File(sys.argv[1],'r')
delta_nu = float(sys.argv[2])

Nx = a['Nx'][0]
Ny = a['Ny'][0]
Nz = a['Nz'][0]
Lx = a['Lx'][0]
Ly = a['Ly'][0]
Lz = a['Lz'][0]
y  = a['collocation_points_y'][()]

print("File:             %s" % sys.argv[1])
print("delta_nu:         %s" % delta_nu)
print("Nx:               %s" % Nx)
print("Ny:               %s" % Ny)
print("Nz:               %s" % Nz)

print("\\Delta{}x^{+}:    %s" % (Lx / Nx / delta_nu))
print("\\Delta{}z^{+}:    %s" % (Lz / Nz / delta_nu))

y /= delta_nu
t = np.abs(y - 10).argmin()
print("y_1^{+}:          %s" % (y[ 1] - y[0]))
print("y_{10}^{+}:       %s" % (y[10] - y[0]))
print("y inside 10^{+}:  %s" % t)
print("y_{%d}^{+}:       %s" % (t-1, y[t-1] - y[0]))
print("y_{%d}^{+}:       %s" % (t, y[t] - y[0]))
print("y_{%d}^{+}:       %s" % (t+1, y[t+1] - y[0]))
