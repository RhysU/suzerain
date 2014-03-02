#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# Shorthand for binary under test for desired operator without statistics
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default
perfect="prun ../perfect_advance ${DECOMP:-} $OPER --statistics_dt=0 --statistics_nt=0"

# These datasets are related to implicit forcing and only are meaningful when
# using --advance_nt=N for N > 1.  They must be ignored for --advance_nt=0.
exclude_datasets="--exclude-path=/bar_f            \
                  --exclude-path=/bar_qb           \
                  --exclude-path=/bar_f_dot_u      \
                  --exclude-path=/bar_Crho         \
                  --exclude-path=/bar_Crhou        \
                  --exclude-path=/bar_CrhoE        \
                  --exclude-path=/bar_Crhou_dot_u"
exclude_datasets=$(echo $exclude_datasets | tr -d '\n' | tr -s ' ')

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

banner "Idempotence of restarting without time advance${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=wisdom.init" # Prepared by test_setup.sh
    $perfect mms0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1
    differ $exclude_datasets mms0.h5 a0.h5
)

banner "Equivalence of a field advanced both with and without a restart${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $perfect mms0.h5 --restart_destination "a#.h5" --advance_nt=2 $WIZ $P \
                     --restart_retain=1
    $perfect a0.h5   --restart_destination "b#.h5" --advance_nt=2 $WIZ $P \
                     --restart_retain=1
    $perfect mms0.h5 --restart_destination "c#.h5" --advance_nt=4 $WIZ $P \
                     --restart_retain=1

    # Ensure simulation time "/t" matches before bothering with anything else
    differ --use-system-epsilon --nan b0.h5 c0.h5 /t
    differ --delta=6e-15 --nan --exclude-path=/t b0.h5 c0.h5
)

banner "Upsample/downsample both homogeneous directions${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $perfect mms0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1                                   \
                     --Nx=$((2*$Nx)) --Nz=$((3*$Nz))
    $perfect a0.h5   --restart_destination "b#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1                                   \
                     --Nx=$((  $Nx)) --Nz=$((  $Nz))
    differ $exclude_datasets mms0.h5 b0.h5
)

banner "Upsample/downsample inhomogeneous direction order${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $perfect mms0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1                                   \
                     --k=$(($k+1))
    $perfect a0.h5   --restart_destination "b#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1                                   \
                     --k=$(($k  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ $exclude_datasets --delta=5e-5 mms0.h5 b0.h5 /rho
    differ $exclude_datasets --delta=3e-4 mms0.h5 b0.h5 /rho_u
    differ $exclude_datasets --delta=5e-5 mms0.h5 b0.h5 /rho_v
    differ $exclude_datasets --delta=6e-5 mms0.h5 b0.h5 /rho_w
    differ $exclude_datasets --delta=7e-4 mms0.h5 b0.h5 /rho_E
)

banner "Upsample/downsample inhomogeneous direction NDOF and htdelta${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $perfect mms0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1                                   \
                     --Ny=$((2*$Ny)) --htdelta=$(($htdelta+1))
    $perfect a0.h5   --restart_destination "b#.h5" --advance_nt=0 $WIZ $P \
                     --restart_retain=1                                   \
                     --Ny=$((  $Ny)) --htdelta=$(($htdelta  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ $exclude_datasets --delta=6e-6 mms0.h5 b0.h5 /rho
    differ $exclude_datasets --delta=1e-4 mms0.h5 b0.h5 /rho_u
    differ $exclude_datasets --delta=7e-6 mms0.h5 b0.h5 /rho_v
    differ $exclude_datasets --delta=3e-5 mms0.h5 b0.h5 /rho_w
    differ $exclude_datasets --delta=2e-4 mms0.h5 b0.h5 /rho_E
)

done
