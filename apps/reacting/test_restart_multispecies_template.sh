#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup_multispecies.sh"

# Shorthand for binary under test for desired operator without statistics
: ${OPER:=} # E.g. '--explicit' or '--implicit' or unset to use default
reacting="prun ../reacting_advance -v $OPER --filter=viscous \
          --statistics_dt=0 --statistics_nt=0"

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
    $reacting multi0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ $P \
                        --restart_retain=1
    differ $exclude_datasets multi0.h5 a0.h5
)

banner "Equivalence of a field advanced both with and without a restart${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $reacting multi0.h5 --restart_destination "a#.h5" --advance_nt=2 $WIZ $P \
                        --restart_retain=1
    $reacting a0.h5     --restart_destination "b#.h5" --advance_nt=2 $WIZ $P \
                        --restart_retain=1
    $reacting multi0.h5 --restart_destination "c#.h5" --advance_nt=4 $WIZ $P \
                        --restart_retain=1

    # Ensure simulation time "/t" matches before bothering with anything else
    differ --use-system-epsilon --nan b0.h5 c0.h5 /t
    differ --delta=6e-15 --nan --exclude-path=/t b0.h5 c0.h5
)

banner "Upsample/downsample both homogeneous directions${OPER:+ ($OPER)}"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $reacting multi0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ $P \
                        --restart_retain=1 --Nx=$((2*$Nx)) --Nz=$((3*$Nz))
    $reacting a0.h5     --restart_destination "b#.h5" --advance_nt=0 $WIZ $P \
                        --restart_retain=1 --Nx=$((  $Nx)) --Nz=$((  $Nz))
    differ $exclude_datasets multi0.h5 b0.h5
)

done
