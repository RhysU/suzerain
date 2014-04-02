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

for dataset in /bar_u_u /bar_u_u_u /bar_T_u /helm /metadata_generated                       \
               /bar_om /bar_om_om /bar_rho_om /bar_rho_om_om /bar_rho2 /bar_p2 /bar_T2      \
               /sym_rho2_grad_u /rho2_grad_T /rho2_u_u /rho2_u_u_u                          \
               /rho2_T_u /rho2_mu /rho2_om /rho2_om_om                                      \
               /maxx_rho  /maxx_rho_E  /maxx_rho_u  /maxx_rho_v  /maxx_rho_w                \
               /maxz_rho  /maxz_rho_E  /maxz_rho_u  /maxz_rho_v  /maxz_rho_w                \
               /minx_rho  /minx_rho_E  /minx_rho_u  /minx_rho_v  /minx_rho_w                \
               /minz_rho  /minz_rho_E  /minz_rho_u  /minz_rho_v  /minz_rho_w                \
               /bar_C2rho /bar_C2rhoE /bar_C2rhou /bar_C2rhou /bar_C2rhow /bar_C2rhou_dot_u
do
    excludes+=" --exclude-path $dataset"
done

"$perfect_initial" initial.h5 --k=6 --Nx=1 --Ny=36 --Nz=1 --htdelta=-3 --Re=3000 --Ma=1.5 --Pr=0.7 --Ly=1 --Re_x=3000 --largo_formulation=temporal --largo_grdelta=1e-1
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/initial.h5" initial.h5
"$perfect_advance" --explicit $(readlink -f initial.h5) --status_nt=1 --advance_dt=0.000919
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/restart00000.h5" restart00000.h5

echo Success for $case
