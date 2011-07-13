#!/bin/bash
set -eu
source "`dirname $0`/test_channel_setup.sh"  # Test infrastructure

for METACASE in NP=1 NP=2; do
eval "$METACASE"

banner "Idempotence of restarting without time advancement"
(
    cd $tmpdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0
    differ mms0.h5 a0.h5
)

banner "Equivalence of a field both with and without a restart"
(
    cd $tmpdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=1
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=1
    prunq ../channel_explicit mms0.h5 --desttemplate "c#.h5" --advance_nt=2
    differ --use-system-epsilon b0.h5 c0.h5
)

banner "Upsample/downsample both homogeneous directions"
(
    cd $tmpdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                      --Nx=$((2*$Nx)) --Nz=$((3*$Nz))
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0 \
                                      --Nx=$((  $Nx)) --Nz=$((  $Nz))
    differ mms0.h5 b0.h5
)

banner "Upsample/downsample inhomogeneous direction order"
(
    cd $tmpdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                      --k=$(($k+1))
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0 \
                                      --k=$(($k  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ --delta=5e-5 mms0.h5 b0.h5 /rho
    differ --delta=3e-4 mms0.h5 b0.h5 /rhou
    differ --delta=5e-5 mms0.h5 b0.h5 /rhov
    differ --delta=6e-5 mms0.h5 b0.h5 /rhow
    differ --delta=7e-4 mms0.h5 b0.h5 /rhoe
)

banner "Upsample/downsample inhomogeneous direction NDOF and htdelta"
(
    cd $tmpdir
    runq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0     \
                                     --Ny=$((2*$Ny)) --htdelta=$(($htdelta+1))
    runq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0     \
                                     --Ny=$((  $Ny)) --htdelta=$(($htdelta  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ --delta=6e-6 mms0.h5 b0.h5 /rho
    differ --delta=1e-4 mms0.h5 b0.h5 /rhou
    differ --delta=7e-6 mms0.h5 b0.h5 /rhov
    differ --delta=3e-5 mms0.h5 b0.h5 /rhow
    differ --delta=2e-4 mms0.h5 b0.h5 /rhoe
)

done
