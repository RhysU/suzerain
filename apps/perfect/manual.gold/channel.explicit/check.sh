#!/bin/bash
set -eu

## A testing script to see if explicit channel advance matches gold solutions.
## Small mismatch is acceptable.  Notice --use-system-epsilon can be
## provided as a command line option so h5diff is less picky.

SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
rmmkcd() { rm -rf "$1" && mkdir -vp "$1" && cd "$1"; }
perfect_init=$(readlink -f perfect_init)
perfect_advance=$(readlink -f perfect_advance)
neno="$SCRIPTDIR/../neno"

case=$(basename "$SCRIPTDIR")
rmmkcd "gold/$case"

"$neno" "$perfect_init" --Ny=32 initial.h5 --htdelta=2 \
    > init.log
"$neno" h5diff -r -c -v --exclude-path /metadata_generated "$@" "$SCRIPTDIR/initial.h5" initial.h5   \
    > initial.diff

"$neno" "$perfect_advance" --explicit initial.h5  --status_nt=1 --advance_dt=0.000919 \
    > advance.log
"$neno" h5diff -r -c -v --exclude-path /metadata_generated "$@" "$SCRIPTDIR/restart0.h5" restart0.h5 \
    > restart0.diff

echo Success for $case
