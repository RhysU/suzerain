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
    --Ny=64
    --htdelta=-1/2   # Flat plate case
    --Ly=1/2         # Shortens wall-normal traversal times
    --bulk_rho=1     # Constant density to any power is one
    --npower=0       # Provides linear velocity profile...
    --bulk_rho_u=0   # ...in conjunction with this setting
    --lower_T=1
    --upper_T=1
    --lower_u=1
    --upper_u=1
    --Ma=1.5         # From upper_u/sqrt(upper_T)
    --Re=1000        # Larger reduces undesirable viscous impact
                     # Try --Re=inf for an Eulerian good time
)
perfect_advance=(
    $(readlink -f perfect_advance)
    initial.h5
    -v
    --explicit
#   --undriven=all         # FIXME To increase coverage
    --advance_dt=5         # TODO  Increase above 1
    --status_dt=0.005
    --statistics_dt=0.05
)
perfect_mean=(
    $(readlink -f perfect_mean)
)

# Common case execution logic
run_case() {
    ${perfect_advance[*]}
}

# Common case post-processing logic
run_postproc() {
    grep bc.upper log.dat > upper.dat
    grep bc.lower log.dat > lower.dat
    shopt -s nullglob
    ${perfect_mean[*]} -o summary.h5 sample*.h5 initial.h5 restart*.h5
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
    ) || echo 'v > 0' >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/quiet_pos_vw
    (
        ${perfect_init[*]} --lower_v=0.01 --lower_w=0.005 --upper_v=0.01 --upper_w=0.005
        run_case
    ) || echo 'v, w > 0' >> $FAILURES
    run_postproc
) &


echo 'Quiet base flow suite two: inflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/quiet_neg_v
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=0 --upper_v=-0.01 --upper_w=0
        run_case
    ) || echo 'v < 0' >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/quiet_neg_vw
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=-0.005 --upper_v=-0.01 --upper_w=-0.005
        run_case
    ) || echo 'v, w < 0' >> $FAILURES
    run_postproc
) &

echo 'Pulse: outflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/acoustic_pos_v
    (
        ${perfect_init[*]} --lower_v=0.01 --lower_w=0 --upper_v=0.01 --upper_w=0 \
                           --acoustic_strength=0.01
        run_case
    ) || echo 'v > 0' >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/entropy_pos_v
    (
        ${perfect_init[*]} --lower_v=0.01 --lower_w=0 --upper_v=0.01 --upper_w=0 \
                           --entropy_strength=0.01
        run_case
    ) || echo 'v > 0' >> $FAILURES
    run_postproc
) &

echo 'Pulse: inflow upper boundary'
echo '######################################################'
(
    rmmkcd ticket2399/acoustic_neg_v
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=0 --upper_v=-0.01 --upper_w=0 \
                           --acoustic_strength=0.01
        run_case
    ) || echo 'v < 0' >> $FAILURES
    run_postproc
) &

(
    rmmkcd ticket2399/entropy_neg_v
    (
        ${perfect_init[*]} --lower_v=-0.01 --lower_w=0 --upper_v=-0.01 --upper_w=0 \
                           --entropy_strength=0.01
        run_case
    ) || echo 'v < 0' >> $FAILURES
    run_postproc
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
