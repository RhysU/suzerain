#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# Shorthand for binary under test for desired operator without statistics
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default
perfect="prun ../perfect_advance $OPER --statistics_dt=0 --statistics_nt=0"

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

# We want to share wisdom across test cases as much as possible to tamp down
# rounding-related discrepancies due to FFT kernel differences (ticket #2515).
WIZ="--plan_wisdom=$(mktemp "$testdir/wisdom.XXXXXX")"

# Prepare pich0.h5 in serial, then in parallel, and ensure both match
# pich0.h5 restart file is used in the tests that follow
banner "Preparation of physical-space version of wave-based test field${OPER:+ ($OPER)}"
(
    cd $testdir
    run ../perfect_advance ich0.h5 ${DECOMP:-} $WIZ $OPER    \
                           --restart_destination "pich#.h5"  \
                           --restart_retain=1                \
                           --advance_nt=0 --restart_physical
    $perfect ich0.h5 --restart_destination "a#.h5" --advance_nt=0 \
                     --restart_retain=1                           \
                     $WIZ ${DECOMP:-} --restart_physical
    differ --delta=5e-16 pich0.h5 a0.h5
)

banner "Idempotence of restarting from physical space without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    $perfect pich0.h5 --restart_destination "a#.h5" --advance_nt=0 \
                      --restart_retain=1                           \
                      ${DECOMP:-} $WIZ $P --restart_physical
    differ --delta=1e-15 --exclude-path /twopoint_kx --exclude-path /twopoint_kz pich0.h5 a0.h5
)

banner "Conversion from physical- to wave-based restart without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    $perfect pich0.h5 --restart_destination "a#.h5"   \
                      --restart_retain=1              \
                      --advance_nt=0 $WIZ ${DECOMP:-}
    differ --delta=5e-15 ich0.h5 a0.h5
)

banner "Equivalence of a field advanced both with and without a physical space restart${OPER:+ ($OPER)}"
# --max_dt option avoids false positives from large implicit
# timesteps necessarily magnifying O(epsilon) restart errors
(
    cd $testdir
    $perfect pich0.h5 --restart_destination "a#.h5" --advance_nt=2 $WIZ $P \
                      --restart_retain=1                                   \
                      ${DECOMP:-} --restart_physical --max_dt=1e-5
    $perfect a0.h5    --restart_destination "b#.h5" --advance_nt=2 $WIZ $P \
                      --restart_retain=1                                   \
                      ${DECOMP:-} --restart_physical --max_dt=1e-5
    $perfect pich0.h5 --restart_destination "c#.h5" --advance_nt=4 $WIZ $P \
                      --restart_retain=1                                   \
                      ${DECOMP:-} --restart_physical --max_dt=1e-5
    differ $exclude_datasets_bar --delta=6e-13 b0.h5 c0.h5
    # Paths like /bar_foo not checked as part of this test
)

done
