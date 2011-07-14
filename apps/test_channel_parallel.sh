#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_channel_setup.sh"

: ${ADVANCE:=--advance_nt=10}

banner "Generating serial result for comparison purposes"
(
    cd $testdir
    runq ../channel_explicit mms0.h5 --desttemplate "serial#.h5" $ADVANCE
)

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=3;P=--Pa=3' 'NP=3;P=--Pb=3' 'NP=4'}; do
NP=
P=
eval "$METACASE"

banner "Equivalence of serial and parallel execution"
(
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" $ADVANCE $P
    differ --use-system-epsilon serial0.h5 a0.h5
)

done
