#!/bin/bash
set -eu

## A testing script to see if implicit=rhome_xyz channel advance matches gold solutions.
## Small mismatch is acceptable.  Notice --use-system-epsilon can be
## provided as a command line option so h5diff is less picky.

SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
rmmkcd() { rm -rf "$1" && mkdir -vp "$1" && cd "$1"; }
perfect_initial=$(readlink -f perfect_initial)
perfect_advance=$(readlink -f perfect_advance)

case=$(basename "$SCRIPTDIR")
rmmkcd "gold/$case"
exec 1> >(tee ./output) 2>&1

#Exclusion of bar_ke post-r43126 represents laziness in updating gold files
#Exclusion of bar_p  post-r43655 represents laziness in updating gold files
excludes="--exclude-path /metadata_generated --exclude-path /bar_ke --exclude-path /bar_p"

"$perfect_initial" --Ny=32 initial.h5 --htdelta=2
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/initial.h5" initial.h5
"$perfect_advance" --implicit=rhome_xyz $(readlink -f initial.h5) --status_nt=1 --advance_dt=0.000919
h5diff -r -c -v $excludes "$@" "$SCRIPTDIR/restart00000.h5" restart00000.h5

echo Success for $case
