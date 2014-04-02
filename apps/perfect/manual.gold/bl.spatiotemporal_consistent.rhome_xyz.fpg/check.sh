#!/bin/bash
set -eu

## A testing script to see if explicit boundary layer advance matches gold solutions.
## Small mismatch is acceptable.  Notice --use-system-epsilon can be
## provided as a command line option so h5diff is less picky.

SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
rmmkcd() { rm -rf "$1" && mkdir -vp "$1" && cd "$1"; }
perfect_initial=$(readlink -f perfect_initial)
perfect_advance=$(readlink -f perfect_advance)

case=$(basename "$SCRIPTDIR")
rmmkcd "gold/$case"
exec 1> >(tee ./output) 2>&1

for dataset in /bar_u_u /bar_u_u_u /bar_T_u /helm /metadata_generated \
               /maxx_rho  /maxx_rho_E  /maxx_rho_u  /maxx_rho_v  /maxx_rho_w \
               /maxz_rho  /maxz_rho_E  /maxz_rho_u  /maxz_rho_v  /maxz_rho_w \
               /minx_rho  /minx_rho_E  /minx_rho_u  /minx_rho_v  /minx_rho_w \
               /minz_rho  /minz_rho_E  /minz_rho_u  /minz_rho_v  /minz_rho_w \
               /bar_C2rho /bar_C2rhoE /bar_C2rhou /bar_C2rhou /bar_C2rhow /bar_C2rhou_dot_u
do
    excludes+=" --exclude-path $dataset"
done

# Initialize the appropriate boundary layer grid lacking the radial flow settings
"$perfect_initial" initial.h5 --k=6 --Nx=1 --Ny=36 --Nz=1 --htdelta=-3 --Re=3000 --Ma=1.5 --Pr=0.7 --Ly=1 --Re_x=3000 --largo_formulation=spatiotemporal_consistent --largo_grdelta=1e-1
h5diff -r -c -v $excludes  "$@" "$SCRIPTDIR/initial.h5" initial.h5

# Overwrite with radial nozzle baseflow conditions 1 meter leeward of laminar CEV stagnation
"$perfect_advance" $(readlink -f initial.h5)          \
                   --restart_destination=initial\#.h5 \
                   --explicit --advance_nt=0          \
                   --cevisslam 1.0
h5diff -r -c -v $excludes  "$@" "$SCRIPTDIR/initial00000.h5" initial00000.h5

# Now advance in time and check the result
"$perfect_advance" --implicit=rhome_xyz $(readlink -f initial00000.h5) --status_nt=1 --advance_dt=0.000919 --max_dt=0.0000919
h5diff -r -c -v $excludes  "$@" "$SCRIPTDIR/restart00000.h5" restart00000.h5

echo Success for $case
