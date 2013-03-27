#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

: ${ADVANCE:=--advance_nt=1 --statistics_dt=0 --statistics_nt=0}
: ${OPER:=--implicit} # E.g. '--explicit' or unset to use default

# We want to share wisdom across test cases as much as possible to tamp down
# rounding-related discrepancies due to FFT kernel differences (ticket #2515)
WIZ="--plan_wisdom=$(mktemp "$testdir/wisdom.XXXXXX")"

banner "Generating default implicit solve result for comparison purposes"
(
    cd $testdir
    ../perfect_advance $OPER $ADVANCE $WIZ mms0.h5 \
                       --restart_destination "a#.h5"
)

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable SPECIFICATION).
for SPECIFICATION in                                     \
    zgbsv                                                \
    zgbsvx,equil=false                                   \
    zgbsvx,equil=true                                    \
    zcgbsvx,reuse=false,aiter=1,siter=-1,diter=5,tolsc=0 \
  # FIXME Ticket #2809 zcgbsvx,reuse=true,aiter=1,siter=-1,diter=5,tolsc=0  \
  # FIXME Ticket #2809 zcgbsvx,reuse=true,aiter=5,siter=25,diter=5,tolsc=0
do

banner "Similarity of solver specification ${SPECIFICATION:+ ($SPECIFICATION)}"
(
    cd $testdir
    ../perfect_advance $OPER $ADVANCE $WIZ --solver=$SPECIFICATION mms0.h5 \
                       --restart_destination "b#.h5"

    # Ensure simulation time "/t" matches before bothering with anything else
    differ --use-system-epsilon --nan a0.h5 b0.h5 /t

    # Only worry about differences in state, not other derived quantities
    differ --delta=5e-15 --nan a0.h5 b0.h5 /rho
    differ --delta=5e-15 --nan a0.h5 b0.h5 /rho_u
    differ --delta=5e-15 --nan a0.h5 b0.h5 /rho_v
    differ --delta=5e-15 --nan a0.h5 b0.h5 /rho_w
    differ --delta=5e-15 --nan a0.h5 b0.h5 /rho_E
)

done
