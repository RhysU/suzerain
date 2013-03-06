#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

: ${ADVANCE:=--advance_nt=10 --status_nt=5 --fluct_percent=10 --fluct_seed=45678}
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default

banner "Generating serial result for comparison purposes${OPER:+ ($OPER)}"
(
    cd $testdir
    run ../reacting_advance $OPER mms0.h5 --restart_destination "serial#.h5" $ADVANCE
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
    prun ../reacting_advance $OPER mms0.h5 --restart_destination "a#.h5" $ADVANCE $P
    # Stricter tolerance performed first for non-/bar_foo quantities
    differ $exclude_datasets_bar --delta=3e-13 --nan serial0.h5 a0.h5
    for dset in $datasets_bar; do
        differ --delta=5e-12 serial0.h5 a0.h5 $dset
    done
)

done
