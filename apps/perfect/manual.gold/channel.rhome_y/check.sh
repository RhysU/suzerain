#!/bin/bash
set -eu

## A testing script to see if implicit=rhome_y channel advance matches gold solutions.
## Small mismatch is acceptable.  Notice --use-system-epsilon can be
## provided as a command line option so h5diff is less picky.

SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
rmmkcd() { rm -rf "$1" && mkdir -vp "$1" && cd "$1"; }
perfect_initial=$(readlink -f perfect_initial)
perfect_advance=$(readlink -f perfect_advance)

case=$(basename "$SCRIPTDIR")
rmmkcd "gold/$case"
exec 1> >(tee ./output) 2>&1

for dataset in /bar_u_u /bar_u_u_u /bar_T_u /bar_ke /bar_p /helm /metadata_generated           \
               /bar_om /bar_om_om /bar_rho_om /bar_rho_om_om                                   \
               /fneg_rho  /fneg_rho_E  /fneg_rho_u  /fneg_rho_v  /fneg_rho_w                   \
               /max_rho   /max_rho_E   /max_rho_u   /max_rho_v   /max_rho_w                    \
               /maxx_rho  /maxx_rho_E  /maxx_rho_u  /maxx_rho_v  /maxx_rho_w                   \
               /maxz_rho  /maxz_rho_E  /maxz_rho_u  /maxz_rho_v  /maxz_rho_w                   \
               /min_rho   /min_rho_E   /min_rho_u   /min_rho_v   /min_rho_w                    \
               /minx_rho  /minx_rho_E  /minx_rho_u  /minx_rho_v  /minx_rho_w                   \
               /minz_rho  /minz_rho_E  /minz_rho_u  /minz_rho_v  /minz_rho_w                   \
               /bar_C2rho /bar_C2rhoE  /bar_C2rhou  /bar_C2rhou  /bar_C2rhow /bar_C2rhou_dot_u
do
    excludes+=" --exclude-path $dataset"
done

"$perfect_initial" --Ny=32 initial.h5 --htdelta=2
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/initial.h5" initial.h5
"$perfect_advance" --implicit=rhome_y $(readlink -f initial.h5) --status_nt=1 --advance_dt=0.000919
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/restart00000.h5" restart00000.h5

echo Success for $case
