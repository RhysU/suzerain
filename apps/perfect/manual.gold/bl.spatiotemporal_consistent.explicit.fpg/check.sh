#!/bin/bash
set -eu

## A testing script to see if explicit boundary layer advance matches gold solutions.
## Small mismatch is acceptable.  Notice --use-system-epsilon can be
## provided as a command line option so h5diff is less picky.

SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
rmmkcd() { rm -rf "$1" && mkdir -vp "$1" && cd "$1"; }
perfect_init=$(readlink -f perfect_init)
perfect_advance=$(readlink -f perfect_advance)

case=$(basename "$SCRIPTDIR")
rmmkcd "gold/$case"
exec 1> >(tee ./output) 2>&1

#Exclusion of bar_ke post-r43126 represents laziness in updating gold files
excludes="--exclude-path /metadata_generated --exclude-path /bar_ke"

# Initialize the appropriate boundary layer grid lacking the radial nozzle settings
"$perfect_init" initial.h5 --k=6 --Nx=1 --Ny=36 --Nz=1 --htdelta=-3 --Re=3000 --Ma=1.5 --Pr=0.7 --Ly=1 --Re_x=3000 --largo_formulation=spatiotemporal_consistent --largo_grdelta=1e-1
h5diff -r -c -v $excludes  "$@" "$SCRIPTDIR/initial.h5" initial.h5

# Overwrite with radial nozzle baseflow conditions about 1 meter leeward of CEV stagnation
"$perfect_advance" $(readlink -f initial.h5)          \
                   --restart_destination=initial\#.h5 \
                   --explicit --advance_nt=0          \
                   --radialflow_Ma0=3.94774352598265  \
                   --radialflow_gam0=1.40926817855057 \
                   --radialflow_rho1=1                \
                   --radialflow_u1=-0.154983103461198 \
                   --radialflow_R1=46.1242981386578
h5diff -r -c -v $excludes  "$@" "$SCRIPTDIR/initial0.h5" initial0.h5

# Now advance in time and check the result
"$perfect_advance" --explicit $(readlink -f initial0.h5) --status_nt=1 --advance_dt=0.000919
h5diff -r -c -v $excludes  "$@" "$SCRIPTDIR/restart0.h5" restart0.h5

echo Success for $case
