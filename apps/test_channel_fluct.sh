#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_channel_setup.sh"

: ${ADVANCE:=--advance_nt=0 --fluctpercent=10 --fluctseed=45678}

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

banner "Reproducibility of adding fluctuations to an existing field"
(
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" $ADVANCE $P
    prunq ../channel_explicit mms0.h5 --desttemplate "b#.h5" $ADVANCE $P
    differ --delta=2e-14 --nan a0.h5 b0.h5
)

done
