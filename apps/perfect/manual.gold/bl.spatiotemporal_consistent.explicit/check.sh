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

for dataset in /twopoint_kx /twopoint_kz /helm /metadata_generated            \
               /fneg_rho  /fneg_rho_E  /fneg_rho_u  /fneg_rho_v  /fneg_rho_w  \
               /max_rho   /max_rho_E   /max_rho_u   /max_rho_v   /max_rho_w   \
               /maxx_rho  /maxx_rho_E  /maxx_rho_u  /maxx_rho_v  /maxx_rho_w  \
               /maxz_rho  /maxz_rho_E  /maxz_rho_u  /maxz_rho_v  /maxz_rho_w  \
               /min_rho   /min_rho_E   /min_rho_u   /min_rho_v   /min_rho_w   \
               /minx_rho  /minx_rho_E  /minx_rho_u  /minx_rho_v  /minx_rho_w  \
               /minz_rho  /minz_rho_E  /minz_rho_u  /minz_rho_v  /minz_rho_w  \
               /bar_Ma /bar_Ma2 /bar_a_u                                      \
               /bar_h02 /bar_H02 /bar_mu2 /bar_nu2 /bar_rho_a /bar_rho2_a     \
               /bar_du_du /bar_rho_du_du /bar_rho2_du_du                      \
               /bar_rho_du /bar_rho2_du                                       \
               /bar_rho_E_u /bar_rho2_E_u /bar_rho_E_E /bar_rho2_E_E          \
               /bar_SrhoE  /bar_Srhou  /bar_Srho  /bar_Srhou_dot_u            \
               /bar_S2rhoE /bar_S2rhou /bar_S2rho /bar_S2rhou_dot_u
do
    excludes+=" --exclude-path $dataset"
done

"$perfect_initial" initial.h5 --k=6 --Nx=1 --Ny=36 --Nz=1 --htdelta=-3 --Re=3000 --Ma=1.5 --Pr=0.7 --Ly=1 --Re_x=3000 --largo_formulation=spatiotemporal_consistent --largo_grdelta=1e-1
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/initial.h5" initial.h5
"$perfect_advance" --explicit $(readlink -f initial.h5) --status_nt=1 --advance_dt=0.000919
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/restart00000.h5" restart00000.h5

echo Success for $case
