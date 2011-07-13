#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_channel_setup.sh"

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=2;P=--Pa=2' 'NP=2;P=--Pb=2'}; do
NP=
P=
eval "$METACASE"

banner "Idempotence of restarting without time advancement"
(
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 $P
    differ mms0.h5 a0.h5
)

banner "Equivalence of a field both with and without a restart"
(
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=1 $P
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=1 $P
    prunq ../channel_explicit mms0.h5 --desttemplate "c#.h5" --advance_nt=2 $P
    differ --use-system-epsilon b0.h5 c0.h5
)

banner "Upsample/downsample both homogeneous directions"
(
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 $P \
                                      --Nx=$((2*$Nx)) --Nz=$((3*$Nz))
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0 $P \
                                      --Nx=$((  $Nx)) --Nz=$((  $Nz))
    differ mms0.h5 b0.h5
)

banner "Upsample/downsample inhomogeneous direction order"
(
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 $P \
                                      --k=$(($k+1))
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0 $P \
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
    cd $testdir
    prunq ../channel_explicit mms0.h5 --desttemplate "a#.h5" --advance_nt=0 $P \
                                     --Ny=$((2*$Ny)) --htdelta=$(($htdelta+1))
    prunq ../channel_explicit a0.h5   --desttemplate "b#.h5" --advance_nt=0 $P \
                                     --Ny=$((  $Ny)) --htdelta=$(($htdelta  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ --delta=6e-6 mms0.h5 b0.h5 /rho
    differ --delta=1e-4 mms0.h5 b0.h5 /rhou
    differ --delta=7e-6 mms0.h5 b0.h5 /rhov
    differ --delta=3e-5 mms0.h5 b0.h5 /rhow
    differ --delta=2e-4 mms0.h5 b0.h5 /rhoe
)

done
