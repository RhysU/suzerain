#!/bin/bash
set -eu

## A testing script to run every gold solution case in parallel.
## Failing test case names are displayed after all have finished.
## Any flags, e.g. --use-system-epsilon, will be propagated.

# Used to build and report failed cases
tempfile() { tempprefix=$(basename "$0"); mktemp /tmp/${tempprefix}.XXXXXX; }
FAILURES=$(tempfile)
trap 'rm -f $FAILURES' EXIT

# Scripts to test are every subdirectory with a check.sh
SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
for script in "$SCRIPTDIR"/*/check.sh
do
    ( "$SCRIPTDIR/neno" "$script" "$@" || echo "$script" "$@" >> $FAILURES ) &
done

# Wait for all subshells to finish
# Output any failures
wait
if test -s $FAILURES
then
    echo 'Failed cases were as follows: '
    nl $FAILURES
fi

# Return success whenever $FAILURES is empty
! test -s $FAILURES
