#!/bin/bash
set -eu

## A testing script to run every gold solution case in parallel.
## Failing test case names are displayed after all have finished.
## Any flags, e.g. --use-system-epsilon, will be propagated.

# Used to build and report failed cases
tempfile() { tempprefix=$(basename "$0"); mktemp "/tmp/${tempprefix}.XXXXXX"; }
FAILURES=$(tempfile)
SUCCESSES=$(tempfile)
trap 'rm -f "$SUCCESSES" "$FAILURES"' EXIT

# Default arguments to pass to scripts whenever non provided by user
defaultargs=--use-system-epsilon

# Scripts to test are every subdirectory with a check.sh
SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
for script in "$SCRIPTDIR"/*/check.sh
do
    ( "$SCRIPTDIR/neno" "$script" "${@-$defaultargs}"        \
      && echo "$script" "${@-defaultargs}" >> "$SUCCESSES"   \
      || echo "$script" "${@-defaultargs}" >> "$FAILURES"  ) &
done

# Wait for all tests to finish
wait

# Output successes and failures
if test $# -lt 1; then
    echo "The following default options were provided to each case:"
else
    echo "The following options were provided to each case:"
fi
declare -i i=1
for arg in "${@-$defaultargs}"; do
    printf $'%6d  %s\n' "$i" "$arg"
    i+=1
done
if test -s "$SUCCESSES"; then
    echo 'Successful cases were as follows: '
    sort "$SUCCESSES" | nl
fi
if test -s "$FAILURES"; then
    echo 'Failed cases were as follows: '
    sort "$FAILURES" | nl
fi

# Return success whenever $FAILURES is empty
! test -s "$FAILURES"
