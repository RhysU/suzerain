#!/bin/bash
set -eu

## A testing script to see if explicit channel advance matches gold solutions.
## Small mismatch is acceptable.  Notice --use-system-epsilon can be
## provided as a command line option so h5diff is less picky.

SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
rmmkcd() { rm -rfv "$1" && mkdir -vp "$1" && cd "$1"; }
perfect_init=$(readlink -f perfect_init)
perfect_advance=$(readlink -f perfect_advance)

case=$(basename "$SCRIPTDIR")
rmmkcd "gold/$case"
"$perfect_init" initial.h5 --k=6 --Nx=1 --Ny=36 --Nz=1 --htdelta=-3 --Re=3000 --Ma=1.5 --Pr=0.7 --Ly=1 --Re_x=3000
"$SCRIPTDIR/../neno" h5diff -r -c -v --exclude-path /metadata_generated "$@" "$SCRIPTDIR/initial.h5" initial.h5   | tee initial.diff
"$perfect_advance" --explicit initial.h5 --status_nt=1 --advance_dt=0.000919
"$SCRIPTDIR/../neno" h5diff -r -c -v --exclude-path /metadata_generated "$@" "$SCRIPTDIR/restart0.h5" restart0.h5 | tee restart0.diff
echo Success for $case
