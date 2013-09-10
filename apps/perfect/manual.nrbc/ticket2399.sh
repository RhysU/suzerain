#!/bin/bash
set -eu

# Shorthand used within the individual cases
rmmkcd() {
    rm -rfv "$1" && mkdir -vp "$1" && cd "$1";
}
perfect_init=(
    $(readlink -f perfect_init)
    initial.h5
    -v
    --restart_physical     # Simplifies any required IC debugging
    --Ny=64
    --htdelta=-1/2         # Flat plate case
    --Ly=1/2               # Shortens wall-normal traversal times
    --bulk_rho=1           # Constant density to any power is one
    --npower=0             # Provides linear velocity profile...
    --bulk_rho_u=0         # ...in conjunction with this setting
    --lower_T=1
    --upper_T=1
    --lower_u=1
    --upper_u=1
    --Ma=1.5               # From upper_u/sqrt(upper_T)
    --Re=1000              # Larger reduces undesirable viscous impact
                           # Try --Re=inf for an Eulerian good time
)
perfect_advance=(
    $(readlink -f perfect_advance)
    initial.h5
    -v
    --explicit             # TODO Eventually use ${OPER:=}
    --undriven=all         # TODO Disabling increases coverage
)
perfect_mean=(
    $(readlink -f perfect_mean)
)

# Common case execution running $1 time units defaulting to some duration
run_case() {
    dt=${1:-2.5}
    ${perfect_advance[*]} "--advance_dt=${dt}"        \
                          "--status_dt=${dt}/1000"    \
                          "--statistics_dt=${dt}/100"
}

# Common case post-processing logic
run_postproc() {

    # Extract 0th, 1st, and 2nd derivatives of conserved state at the boundaries
    grep bc.lower.d0 log.dat > lower.d0
    grep bc.lower.d1 log.dat > lower.d1
    grep bc.lower.d2 log.dat > lower.d2
    grep bc.upper.d0 log.dat > upper.d0
    grep bc.upper.d1 log.dat > upper.d1
    grep bc.upper.d2 log.dat > upper.d2

    # Post-process the entire collection of sample and restart files
    shopt -s nullglob
    ${perfect_mean[*]} -f summary.dat -o summary.h5 \
                       sample*.h5 initial.h5 restart*.h5

    # Display (hopefully) logarithmic decay of signal versus simulation time
    if hash gplot 2>/dev/null; then
        gplot -p L2.mean.png -lc -f i=6:10 -x t -y L2 L2.mean.dat using 4:i with lines
    fi

    # Summarize the overall field behavior in a manner that aids visual debugging
    if hash gplot 2>/dev/null; then
        gplot -3M -p bar_p.png   -t bar_p   -x t -y y summary.dat using '"t":"y":"bar_p"'
        gplot -3M -p bar_rho.png -t bar_rho -x t -y y summary.dat using '"t":"y":"bar_rho"'
        gplot -3M -p bar_T.png   -t bar_T   -x t -y y summary.dat using '"t":"y":"bar_T"'
        gplot -3M -p bar_v.png   -t bar_v   -x t -y y summary.dat using '"t":"y":"bar_v"'
    fi
}

# Used to build and report failed cases
tempfile() { tempprefix=$(basename "$0"); mktemp /tmp/${tempprefix}.XXXXXX; }
FAILURES=$(tempfile)
trap 'rm -f $FAILURES' EXIT

echo 'Quiet base flow suite one: outflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/quiet_pos_v
    (
        ${perfect_init[*]} --lower_v=0.01 --lower_w=0 --upper_v=0.01 --upper_w=0
        run_case
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/quiet_pos_vw
    (
        ${perfect_init[*]} --lower_v=0.01 --lower_w=0.005 --upper_v=0.01 --upper_w=0.005
        run_case
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &


echo 'Quiet base flow suite two: inflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/quiet_neg_v
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=0 --upper_v=-0.01 --upper_w=0
        run_case
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/quiet_neg_vw
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=-0.005 --upper_v=-0.01 --upper_w=-0.005
        run_case
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &

echo 'Pulse: outflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/acoustic_pos_v
    (
        ${perfect_init[*]} --lower_v=0.01 --lower_w=0 --upper_v=0.01 --upper_w=0 \
                           --acoustic_strength=0.1
        run_case 1
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/entropy_pos_v
    (
        ${perfect_init[*]} --lower_v=1.0 --lower_w=0 --upper_v=1.0 --upper_w=0 \
                           --entropy_strength=0.01
        run_case 1
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &

echo 'Pulse: inflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/acoustic_neg_v
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=0 --upper_v=-0.01 --upper_w=0 \
                           --acoustic_strength=0.1
        run_case 1
    ) || basename `pwd` >> $FAILURES
    run_postproc
) &

(
    # rmmkcd ticket2399/entropy_neg_v
    # Entropy waves do not hit the upper boundary for negative v
) &

# Attach a plot of some it-bounces-back-and-forth measurement as a baseline, or some other diagnostic procedure that I'll repeat. Make this scriptable so that I may re-run these cases quickly after any code changes.
# What about for v_inflow > 0, v_outflow > 0?
# What about for v_inflow < 0, v_outflow < 0?
# Does my acoustic pulse in perfect_init do what I expect in 1D with the upper NRBC?
# For constant reference v_inflow > 0, v_outflow > 0?
# For constant reference v_inflow < 0, v_outflow < 0?
# For variable reference v_inflow > 0, v_outflow > 0?
# For variable reference v_inflow < 0, v_outflow < 0?
# Can I retain this behavior while adding driving forcing to the freestream?
# Repeat for the entropy pulse, working through a new copy of whatever this final list becomes.

# Wait for all subshells to finish
wait

# Output any failures
if test -s $FAILURES
then
    echo 'Failed cases were as follows: '
    nl $FAILURES
fi

# Return success whenever $FAILURES is empty
! test -s $FAILURES
