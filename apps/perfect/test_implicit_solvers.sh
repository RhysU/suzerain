#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

: ${ADVANCE:=--advance_nt=1 --statistics_dt=0 --statistics_nt=0}
: ${OPER:=--implicit=rhome_xyz} # Implicit in three directions is best test

# We want to share wisdom across test cases as much as possible to tamp down
# rounding-related discrepancies due to FFT kernel differences (ticket #2515)
WIZ="--plan_wisdom=$(mktemp "$testdir/wisdom.XXXXXX")"

banner "Generating default implicit solve result for comparison purposes"
(
    cd $testdir
    ../perfect_advance $OPER $ADVANCE $WIZ mms0.h5   \
                       --restart_destination "a#.h5" \
                       --restart_retain=1
)

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable SPECIFICATION).
for SPECIFICATION in                                     \
    zgbsv                                                \
    zgbsvx,equil=false                                   \
    zgbsvx,equil=true                                    \
    zcgbsvx,reuse=false,aiter=1,siter=-1,diter=5,tolsc=0 \
    zcgbsvx,reuse=true,aiter=1,siter=-1,diter=5,tolsc=0  \
    zcgbsvx,reuse=true,aiter=5,siter=25,diter=5,tolsc=0
do

banner "Similarity of solver specification ${SPECIFICATION:+ ($SPECIFICATION)}"
(
    cd $testdir
    ../perfect_advance $OPER ${DECOMP:-} $ADVANCE $WIZ --solver=$SPECIFICATION  \
                       mms0.h5 --restart_destination "b#.h5" --restart_retain=1

    # Ensure simulation time "/t" matches before bothering with anything else
    differ --use-system-epsilon --nan a0.h5 b0.h5 /t

    # Only worry about differences in state, not other derived quantities
    differ --delta=6e-15 --nan a0.h5 b0.h5 /rho
    differ --delta=6e-15 --nan a0.h5 b0.h5 /rho_u
    differ --delta=6e-15 --nan a0.h5 b0.h5 /rho_v
    differ --delta=6e-15 --nan a0.h5 b0.h5 /rho_w
    differ --delta=1e-14 --nan a0.h5 b0.h5 /rho_E
)

done
