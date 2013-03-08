#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# Shorthand for binary under test for desired operator without statistics
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default
reacting="prun ../reacting_advance $OPER --statistics_dt=0 --statistics_nt=0"

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

# Prepare pmms0.h5 in serial, then in parallel, and ensure both match
# pmms0.h5 restart file is used in the tests that follow
banner "Preparation of physical-space version of wave-based test field${OPER:+ ($OPER)}"
(
    cd $testdir
    run ../reacting_advance $OPER mms0.h5 --restart_destination "pmms#.h5" \
                                         --advance_nt=0 --restart_physical
    $reacting mms0.h5 --restart_destination "a#.h5" --advance_nt=0 \
                     --restart_physical
    differ --delta=5e-16 pmms0.h5 a0.h5
)

banner "Idempotence of restarting from physical space without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    $reacting pmms0.h5 --restart_destination "a#.h5" --advance_nt=0 $P \
                      --restart_physical
    #differ --delta=1e-15 pmms0.h5 a0.h5
    differ --relative=5e-14 pmms0.h5 a0.h5
)

banner "Conversion from physical- to wave-based restart without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    $reacting pmms0.h5 --restart_destination "a#.h5" --advance_nt=0
    
    # original: absolute tolerance too tight because now working with
    # dimensional variables
    #
    #differ --delta=5e-15 mms0.h5 a0.h5

    # relative tolerance doesn't work well because many of the entries
    # we want to check are small numbers computed from operations on
    # big numbers, so relative errors are large
    #
    #differ --relative=5e-14 mms0.h5 a0.h5

    # so... use a large absolute tolerance for now 
    #
    # FIXME: Find a better solution.  Want to be able to find
    # differences that are larger than some tolerance in *both* a
    # relative and absolute sense, but h5diff doesn't provide this
    # automatically.
    #
    # NOTE: Tolerance based on rho_E * 1e-15
    #
    differ --delta=2e-10 mms0.h5 a0.h5
)

banner "Equivalence of a field advanced both with and without a physical space restart${OPER:+ ($OPER)}"
# --max_dt option avoids false positives from large implicit
# timesteps necessarily magnifying O(epsilon) restart errors
(
    cd $testdir
    $reacting pmms0.h5 --restart_destination "a#.h5" --advance_nt=2 $P \
                      --restart_physical --max_dt=1e-5
    $reacting a0.h5    --restart_destination "b#.h5" --advance_nt=2 $P \
                      --restart_physical --max_dt=1e-5
    $reacting pmms0.h5 --restart_destination "c#.h5" --advance_nt=4 $P \
                      --restart_physical --max_dt=1e-5
    differ $exclude_datasets_bar --delta=6e-13 b0.h5 c0.h5
    # Paths like /bar_foo not checked as part of this test
)

done
