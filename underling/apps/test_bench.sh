#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# Do unto others...
renice +7 -p $$

: ${FIELD:=--field-global=19x23x29}  # Field size for tests
: ${CHECK:=--check=1e-15}            # Check tolerance (scaled by FIELD!)
: ${PATIENCE:=--estimate}            # Planning rigor
: ${TEE=tee -a $testdir/output}      # Captures last command output to file
: ${GREP=grep error}                 # Filter output displayed

# Shorthand
bench="./underling_bench $PATIENCE"

# Run each test case in this file under the following circumstances
# (which can be overridden by providing the environment variable METACASES).
for METACASE in ${METACASES:= 'NP=1' 'NP=3;P=--dims=3x0' 'NP=3;P=--dims=0x3' 'NP=4'}
do

NP=
P=
eval "$METACASE"

banner $LINENO "Basic out-of-place transposes"
(
    for h in 1 2 5; do
        prun $bench $FIELD $P $CHECK \
            --howmany=$h |$TEE|$GREP 2>&1
        test ${PIPESTATUS[0]} -eq 0 || exit 1
    done
)

banner $LINENO "Basic out-of-place transposes with TRANSPOSED_LONG_N{0,2}"
(
    for h in 1 2 5; do
        prun $bench $FIELD $P $CHECK --trans \
            --howmany=$h |$TEE|$GREP 2>&1
        test ${PIPESTATUS[0]} -eq 0 || exit 1
    done
)

banner $LINENO "Basic in-place transposes"
(
    for h in 1 2 5; do
        prun $bench $FIELD $P $CHECK --mpi-in-place \
            --howmany=$h |$TEE|$GREP 2>&1
        test ${PIPESTATUS[0]} -eq 0 || exit 1
    done
)

banner $LINENO "Basic in-place transposes with TRANSPOSED_LONG_N{0,2}"
(
    for h in 1 2 5; do
        prun $bench $FIELD $P $CHECK --mpi-in-place --trans \
            --howmany=$h |$TEE|$GREP 2>&1
        test ${PIPESTATUS[0]} -eq 0 || exit 1
    done
)

for fft in ccc ccr icc icr iic iir
do

    banner $LINENO "Transpose out-of-place and transform out-of-place ($fft)"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose in-place and transform out-of-place ($fft)"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --mpi-in-place --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose out-of-place and transform in-place ($fft)"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --fft-in-place --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose in-place and transform in-place ($fft)"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --mpi-in-place --fft-in-place --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

done

for fft in ccc icc iic
do

    banner $LINENO "Transpose out-of-place and transform out-of-place ($fft) with TRANSPOSED_LONG_N{0,2}"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --trans --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose in-place and transform out-of-place ($fft) with TRANSPOSED_LONG_N{0,2}"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --mpi-in-place --trans --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose out-of-place and transform in-place ($fft) with TRANSPOSED_LONG_N{0,2}"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --fft-in-place --trans --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose in-place and transform in-place ($fft) with TRANSPOSED_LONG_N{0,2}"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --mpi-in-place --fft-in-place --trans --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

done

for fft in ccc ccr icc icr iic iir
do
    banner $LINENO "Transpose out-of-place and packed transform out-of-place ($fft)"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --pack --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose in-place and packed transform out-of-place ($fft)"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --mpi-in-place --pack --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose out-of-place and packed transform out-of-place ($fft) with TRANSPOSED_LONG_N{0,2}"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --pack --trans --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

    banner $LINENO "Transpose in-place and packed transform out-of-place ($fft) with TRANSPOSED_LONG_N{0,2}"
    (
        for h in 2 6; do
            prun $bench $FIELD $P $CHECK --mpi-in-place --pack --trans --$fft \
                --howmany=$h |$TEE|$GREP 2>&1
            test ${PIPESTATUS[0]} -eq 0 || exit 1
        done
    )

done

done
