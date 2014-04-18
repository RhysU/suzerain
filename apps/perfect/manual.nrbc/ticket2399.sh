#!/bin/bash
set -eu

## A testing script for semi-automated investigation of non-reflecting boundary
## conditions.  See discussion within Redmine tickets #2399, #2957, and #2979.

## Somewhat or entirely unresolved details (with commentary):
#    1) Attach a plot of some it-bounces-back-and-forth measurement as a
#       baseline, or some other diagnostic procedure that I'll repeat. Make
#       this scriptable so that I may re-run these cases quickly after changes.
#
#       (Function run_postproc outputs a plot of L2 for all variables.)
#       (This script is the repeatable mechanism.)
#
#    2) What about for v_inflow > 0, v_outflow > 0?
#    3) What about for v_inflow < 0, v_outflow < 0?
#    4) Does my acoustic pulse in perfect_initial do what I expect in 1D with
#       the upper NRBC?
#
#       (Yes, yes, and yes.)
#
#    5) For constant reference v_inflow > 0, v_outflow > 0?
#    6) For constant reference v_inflow < 0, v_outflow < 0?
#
#       (Yes and yes.)
#
#    7) For variable reference v_inflow > 0, v_outflow > 0?
#    8) For variable reference v_inflow < 0, v_outflow < 0?
#
#       (No and no.  Dropping variable NRBC reference ideas for now.)
#
#    9) Can I retain this behavior while adding driving forcing to the
#       freestream?
#
#       (Seemingly so, but a true test requires more wavenumbers and
#        a pulse riding on the non-zero-zero wavenumber content.
#        Punting for now.)
#
#   10) Repeat for the entropy pulse, working through a new copy of whatever
#       this final list becomes.
#
#       (Entropy pulses have been problematic.  Punting on testing for now.)
#

# Shorthand used within the individual cases
rmmkcd() {
    rm -rfv "$1" && mkdir -vp "$1" && cd "$1";
}
perfect_initial=(
    $(readlink -f perfect_initial)
    initial.h5
    -v
    --restart_physical     # Simplifies any required IC debugging
    --Ny=64
    --htdelta=-1           # Flat plate case
    --Ly=1/4               # Shortens wall-normal traversal times
    --bulk_rho=1           # Constant density to any power is one
    --npower=0             # Provides linear velocity profile...
    --bulk_rho_u=0         # ...in conjunction with this setting
    --lower_T=1
    --upper_T=1
    --lower_u=1
    --upper_u=1
    --Ma=1.5               # From upper_u/sqrt(upper_T)
    --Re=500               # Larger reduces undesirable viscous impact
                           # Try --Re=inf for an Eulerian good time
)
perfect_advance=(
    $(readlink -f perfect_advance)
    initial.h5
    --undriven=all         # Driving forces interfere with wall-normal pulses
    ${OPER:+$OPER}         # Permit runtime selection of time advance
)
perfect_summary=(
    $(readlink -f perfect_summary)
)

# Common case execution running $1 time units defaulting to some duration
run_case() {
    dt=${1:-1}
    ${perfect_advance[*]} "--advance_dt=${dt}"        \
                          "--status_dt=${dt}/1000"    \
                          "--statistics_dt=${dt}/100" \
                          "$@"
}

# Common case post-processing logic
run_postproc() {

    # Plot 0th, 1st, and 2nd derivatives of conserved state at the boundaries
    if hash gplot 2>/dev/null; then
        gplot -o lower.d0.png -g bc.lower.d0 -x t -c -f i=6:10 using 4:i w l ::: bc.dat
        gplot -o lower.d1.png -g bc.lower.d1 -x t -c -f i=6:10 using 4:i w l ::: bc.dat
        gplot -o lower.d2.png -g bc.lower.d2 -x t -c -f i=6:10 using 4:i w l ::: bc.dat
        gplot -o upper.d0.png -g bc.upper.d0 -x t -c -f i=6:10 using 4:i w l ::: bc.dat
        gplot -o upper.d1.png -g bc.upper.d1 -x t -c -f i=6:10 using 4:i w l ::: bc.dat
        gplot -o upper.d2.png -g bc.upper.d2 -x t -c -f i=6:10 using 4:i w l ::: bc.dat
    fi

    # Display (hopefully) logarithmic decay of signal versus simulation time
    if hash gplot 2>/dev/null; then
        gplot -o L2.png -g state.L2 -lc -f i=6:10 -x t -y L2 using 4:i with lines ::: state.dat
    fi

    # Post-process the entire collection of sample and restart files
    shopt -s nullglob
    ${perfect_summary[*]} -f summary.dat -o summary.h5 \
                       sample*.h5 initial.h5 restart*.h5

    # Summarize the overall field behavior in a manner that aids visual debugging
    if hash gplot 2>/dev/null; then
        gplot -3M -o bar_p.png   -t bar_p   -x t -y y using '"t":"y":"bar_p"'   ::: summary.dat
        gplot -3M -o bar_rho.png -t bar_rho -x t -y y using '"t":"y":"bar_rho"' ::: summary.dat
        gplot -3M -o bar_T.png   -t bar_T   -x t -y y using '"t":"y":"bar_T"'   ::: summary.dat
        gplot -3M -o bar_v.png   -t bar_v   -x t -y y using '"t":"y":"bar_v"'   ::: summary.dat
    fi
}

# Used to build and report failed cases
tempfile() { tempprefix=$(basename "$0"); mktemp /tmp/${tempprefix}.XXXXXX; }
FAILURES=$(tempfile)
trap 'rm -f $FAILURES' EXIT

echo 'Quiet base flow suite one: outflow upper boundary'
echo '######################################################'
(
    # Ticket #2399 showed problems with subsonic viscous outflow.
    # Ticket #2983 failed in an attempt to fix them quickly.
    # FIXME A thorough, IMEX-friendly solution requires some thought.
    # This test case will misbehave until the #2983 and #2957 are resolved.
    case="quiet_pos_v"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=0.1 --lower_w=0 \
                              --upper_v=0.1 --upper_w=0
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

(
    case="quiet_pos_vw"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=0.1 --lower_w=0.05 \
                              --upper_v=0.1 --upper_w=0.05
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &


echo 'Quiet base flow suite two: inflow upper boundary'
echo '######################################################'
(
    case="quiet_neg_v"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=-0.1 --lower_w=0 \
                              --upper_v=-0.1 --upper_w=0
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

(
    case="quiet_neg_vw"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=-0.1 --lower_w=-0.05 \
                              --upper_v=-0.1 --upper_w=-0.05
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

echo 'Acoustic pulse: outflow upper boundary'
echo '######################################################'
(
    case="acoustic_pos_v"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=0.01 --lower_w=0 \
                              --upper_v=0.01 --upper_w=0 \
                              --acoustic_strength=0.1
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

(
    case="acoustic_pos_vw"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=0.1 --lower_w=0.05 \
                              --upper_v=0.1 --upper_w=0.05 \
                              --acoustic_strength=0.1
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

echo 'Acoustic pulse: inflow upper boundary'
echo '######################################################'
(
    case="acoustic_neg_v"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=-0.01 --lower_w=0 \
                              --upper_v=-0.01 --upper_w=0 \
                              --acoustic_strength=0.1
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

(
    case="acoustic_neg_vw"
    rmmkcd "ticket2399/$case"
    (
        ${perfect_initial[*]} --lower_v=-0.1 --lower_w=-0.05 \
                              --upper_v=-0.1 --upper_w=-0.05 \
                              --acoustic_strength=0.1
        run_case
    ) || echo "$case" >> $FAILURES
    run_postproc
) &

## TODO See item 10 above
#
# echo 'Entropy pulse: outflow upper boundary'
# echo '######################################################'

## An entropy pulse makes no sense for an inflow upper boundary
## as convection does not drive it towards the boundary in question.
#
# echo 'Entropy pulse: inflow upper boundary'
# echo '######################################################'

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
