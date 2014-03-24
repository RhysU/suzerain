#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# FIXME Ticket 2065 use 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2' in METACASES
# FIXME Ticket 2065 using run/prun and -v -V options to aid debugging

# Shorthand for serial/parallel run truncating to only 1D problem
# Such runs certainly stress the MPI pencil decomposition routines
# TODO: Switch back to implicit once that capability is resurrected
s_reacting="run  ../reacting_advance --explicit --Nx=1 --Nz=1 --restart_retain=1"
p_reacting="prun ../reacting_advance --explicit --Nx=1 --Nz=1 --restart_retain=1"

# These datasets are related to implicit forcing and only are meaningful when
# using --advance_nt=N for N > 1.  They must be ignored for --advance_nt=0.
exclude_datasets="--exclude-path=/bar_f             \
                  --exclude-path=/bar_qb            \
                  --exclude-path=/bar_f_dot_u       \
                  --exclude-path=/bar_Crho          \
                  --exclude-path=/bar_Crhou         \
                  --exclude-path=/bar_CrhoE         \
                  --exclude-path=/bar_Crhou_dot_u   \
                  --exclude-path=/bar_C2rho         \
                  --exclude-path=/bar_C2rhou        \
                  --exclude-path=/bar_C2rhoE        \
                  --exclude-path=/bar_C2rhou_dot_u"
exclude_datasets=$(echo $exclude_datasets | tr -d '\n' | tr -s ' ')

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1'                                }; do
NP=
P=
eval "$METACASE"

banner "Idempotence of serial versus degenerate parallel without time advance"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $s_reacting mms0.h5 --restart_destination "a#.h5" --advance_nt=0 $WIZ
    $p_reacting mms0.h5 --restart_destination "b#.h5" --advance_nt=0 $WIZ $P
    differ --use-system-epsilon $exclude_datasets a0.h5 b0.h5
)

banner "Idempotence of serial versus degenerate parallel with time advance"
(
    cd $testdir
    WIZ="--plan_wisdom=$(mktemp wisdom.XXXXXX)"
    $s_reacting mms0.h5 -v --restart_destination "a#.h5" --advance_nt=3 $WIZ \
                       --status_nt=1 --evmagfactor=0.2
    $p_reacting mms0.h5 -v --restart_destination "b#.h5" --advance_nt=3 $WIZ \
                       --status_nt=1 --evmagfactor=0.2 $P
    differ --use-system-epsilon a0.h5 b0.h5
)

done
