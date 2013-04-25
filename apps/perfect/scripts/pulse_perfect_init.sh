#!/bin/bash

# Script to generate field with an low amplitude traveling pulse 
# to test the nonreflecting boundary condition

declare -ir Nx=1
declare -ir Ny=100
declare -ir k=6
declare -ir htdelta=1
declare -ir Nz=1
bulk_rho="NaN"

filename="pulse_perfect_physical.h5" 
rm ${filename}

# 
perfect_init ${filename}  \
                     --Nx=$Nx --Ny=$Ny --k=$k --htdelta=$htdelta --Nz=$Nz \
                     --restart_physical --bulk_rho=${bulk_rho} 

# Run python script to impose pulse field
python pulse_perfect.py





