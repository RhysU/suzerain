#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

: ${ADVANCE:=--advance_nt=10 --status_nt=5 --fluct_percent=10 --fluct_seed=45678}
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default

# We want to share wisdom across test cases as much as possible to tamp down
# rounding-related discrepancies due to FFT kernel differences (ticket #2515)
WIZ="--plan_wisdom=$(mktemp "$testdir/wisdom.XXXXXX")"

banner "Generating serial result for comparison purposes${OPER:+ ($OPER)}"
(
    cd $testdir
    run ../perfect_advance ich0.h5 $OPER --restart_destination "serial#.h5" \
                                         --restart_retain=1                 \
                                         $ADVANCE $WIZ ${DECOMP:-}
)

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=3;P=--Pa=3' 'NP=3;P=--Pb=3' 'NP=4'}; do
NP=
P=
eval "$METACASE"

banner "Equivalence of serial and parallel execution${OPER:+ ($OPER)}"
(
    cd $testdir
    prun ../perfect_advance ich0.h5 $OPER --restart_destination "a#.h5" \
                                          --restart_retain=1            \
                                          $ADVANCE $P $WIZ ${DECOMP:-}
    # Stricter tolerance performed first for non-/bar_foo quantities
    differ $exclude_datasets_bar \
           $exclude_datasets_maxminloc --delta=3e-13 --nan serial0.h5 a0.h5
    for dset in $datasets_bar; do
        differ --delta=5e-12 serial0.h5 a0.h5 $dset
    done
)

done
