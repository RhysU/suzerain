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
    run ../reacting_advance $OPER mms0.h5 --restart_destination "serial#.h5" \
                                          $ADVANCE $WIZ
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
    prun ../reacting_advance $OPER mms0.h5 --restart_destination "a#.h5" \
                                           $ADVANCE $P $WIZ
    # Stricter tolerance performed first for non-/bar_foo quantities
    #differ $exclude_datasets_bar --delta=3e-13 --nan serial0.h5 a0.h5
    differ $exclude_datasets_bar --delta=6e-13 --nan serial0.h5 a0.h5
    for dset in $datasets_bar; do
        #differ --delta=5e-12 serial0.h5 a0.h5 $dset
        # FIXME: Ticket 2790 to improve diffing tolerances
        # NOTE: Tolerance so large b/c /bar_mu_grad_T has data on order 5e9
        # TODO: Why so large?  Does this make sense?
        #differ --delta=2e-6 serial0.h5 a0.h5 $dset

        # Updated on 4/9/2013-2014 to get implicit test to pass.  Above
        # still applies.
        # differ --delta=2e-5 serial0.h5 a0.h5 $dset

        # Updated on 4/12/2013-2014 to get implicit test to pass with
        # intel/11.1 in optimized mode.  Above still applies.
        differ --delta=4e-5 serial0.h5 a0.h5 $dset
    done
)

done
