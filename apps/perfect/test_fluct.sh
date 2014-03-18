#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

: ${ADVANCE:=--advance_nt=0 --fluct_percent=10 --fluct_seed=45678  \
                            --fluct_kxfrac=0.5: --fluct_kzfrac=:0.5}

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

banner "Reproducibility of adding fluctuations to an existing field"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    prun ../perfect_advance ich0.h5 --restart_destination "a#.h5" \
                                    --restart_retain=1            \
                                    $ADVANCE $WIZ ${DECOMP:-} $P
    prun ../perfect_advance ich0.h5 --restart_destination "b#.h5" \
                                    --restart_retain=1            \
                                    $ADVANCE $WIZ ${DECOMP:-} $P
    differ --delta=2e-14 --nan a0.h5 b0.h5
)

done
