#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_channel_setup.sh"

# Shorthand
explicit="prunq ../channel_explicit"

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

# Prepare pmms0.h5 in serial, then in parallel, and ensure both match
# pmms0.h5 restart file is used in the tests that follow
banner "Preparation of physical-space version of wave-based test field"
(
    cd $testdir
    runq ../channel_explicit mms0.h5 --desttemplate "pmms#.h5" --advance_nt=0 \
                                     --restart_physical
    $explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 --restart_physical
    differ pmms0.h5 a0.h5
)

banner "Idempotence of restarting from physical space without time advance"
(
    cd $testdir
    $explicit pmms0.h5 --desttemplate "a#.h5" --advance_nt=0 $P --restart_physical
    differ --delta=1e-15 --nan pmms0.h5 a0.h5
)

banner "Conversion from physical- to wave-based restart without time advance"
(
    cd $testdir
    $explicit pmms0.h5 --desttemplate "a#.h5" --advance_nt=0
    differ --delta=3e-15 --nan mms0.h5 a0.h5
)

# BORKED?
# banner "Equivalence of a field both with and without a physical space restart"
# (
#     cd $testdir
#     $explicit pmms0.h5 --desttemplate "a#.h5" --advance_nt=1 $P --restart_physical
#     $explicit a0.h5    --desttemplate "b#.h5" --advance_nt=1 $P --restart_physical
#     $explicit pmms0.h5 --desttemplate "c#.h5" --advance_nt=2 $P --restart_physical
#     differ --delta=7e-16 --nan b0.h5 c0.h5
# )

# BORKED?
# banner "Idempotence of restarting without time advance via physical space"
# (
#     cd $testdir
#     $explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 $P --restart_physical
#     $explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0 $P
#     differ --delta=3e-15 mms0.h5 b0.h5
# )

done
