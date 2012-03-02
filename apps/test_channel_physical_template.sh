#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_channel_setup.sh"

# Shorthand for binary under test for desired operator without statistics
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default
channel="prun ../channel $OPER --statistics_dt=0 --statistics_nt=0"

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
    run ../channel $OPER mms0.h5 --restart_destination "pmms#.h5" \
                                 --advance_nt=0 --restart_physical
    $channel mms0.h5 --restart_destination "a#.h5" --advance_nt=0 \
                     --restart_physical
    differ pmms0.h5 a0.h5
)

banner "Idempotence of restarting from physical space without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    $channel pmms0.h5 --restart_destination "a#.h5" --advance_nt=0 $P \
                      --restart_physical
    differ --delta=1e-15 --nan pmms0.h5 a0.h5
)

banner "Conversion from physical- to wave-based restart without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    $channel pmms0.h5 --restart_destination "a#.h5" --advance_nt=0
    differ --delta=3e-15 --nan mms0.h5 a0.h5
)

banner "Equivalence of a field advanced both with and without a physical space restart${OPER:+ ($OPER)}"
# --max_dt option avoids false positives from large implicit
# timesteps necessarily magnifying O(epsilon) restart errors
(
    cd $testdir
    $channel pmms0.h5 --restart_destination "a#.h5" --advance_nt=2 $P \
                      --restart_physical --max_dt=1e-5
    $channel a0.h5    --restart_destination "b#.h5" --advance_nt=2 $P \
                      --restart_physical --max_dt=1e-5
    $channel pmms0.h5 --restart_destination "c#.h5" --advance_nt=4 $P \
                      --restart_physical --max_dt=1e-5
    differ_exclude $exclude_datasets_bar --delta=6e-13 --relative=6e-13 --nan b0.h5 c0.h5
    # Paths like /bar_foo not checked as part of this test
)

done
