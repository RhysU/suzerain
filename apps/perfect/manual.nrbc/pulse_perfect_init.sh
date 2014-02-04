#!/bin/bash
# Script to generate field with an low amplitude traveling pulse
# to test the nonreflecting boundary condition
set -eu

Nx=1
Ny=96
k=6
htdelta=1
Nz=1
bulk_rho=NaN

# Generate a restart file and then modify it with an initial pulse
filename="pulse_perfect_physical.h5"
../perfect_initial --clobber ${filename}                             \
                --Nx=$Nx --Ny=$Ny --k=$k --htdelta=$htdelta --Nz=$Nz \
                --restart_physical --bulk_rho=${bulk_rho}
./pulse_perfect.py "$filename"
