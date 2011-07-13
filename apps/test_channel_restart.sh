#!/bin/bash
set -e

source "${SRCDIR:-.}/test_channel_setup.sh"  # Test infrastructure

banner "Idempotence of restarting without time advancement"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0
    differ initial.h5 a0.h5
)

banner "Equivalence of a field both with and without a restart"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=1
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=1
    runq ../channel_explicit initial.h5 --desttemplate "c#.h5" --advance_nt=2
    differ --use-system-epsilon b0.h5 c0.h5
)

banner "Upsample/downsample both homogeneous directions"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                        --Nx=$((2*$Nx)) --Nz=$((3*$Nz))
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=0 \
                                        --Nx=$((  $Nx)) --Nz=$((  $Nz))
    differ initial.h5 b0.h5
)

banner "Upsample/downsample inhomogeneous direction order"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                        --k=$(($k+1))
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=0 \
                                        --k=$(($k  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ --delta=5e-5 initial.h5 b0.h5 /rho
    differ --delta=3e-4 initial.h5 b0.h5 /rhou
    differ --delta=5e-5 initial.h5 b0.h5 /rhov
    differ --delta=6e-5 initial.h5 b0.h5 /rhow
    differ --delta=7e-4 initial.h5 b0.h5 /rhoe
)

banner "Upsample/downsample inhomogeneous direction NDOF and htdelta"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                        --Ny=$((2*$Ny)) --htdelta=$(($htdelta + 1))
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=0 \
                                        --Ny=$((  $Ny)) --htdelta=$(($htdelta    ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ --delta=6e-6 initial.h5 b0.h5 /rho
    differ --delta=1e-4 initial.h5 b0.h5 /rhou
    differ --delta=7e-6 initial.h5 b0.h5 /rhov
    differ --delta=3e-5 initial.h5 b0.h5 /rhow
    differ --delta=2e-4 initial.h5 b0.h5 /rhoe
)
