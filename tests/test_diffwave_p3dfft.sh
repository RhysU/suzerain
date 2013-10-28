#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Requires AM_TESTS_FD_REDIRECT in Makefile.am
alert() { if test -t 9; then echo WARN: "$@" >&9; else echo WARN: "$@"; fi; }

if ! [ -x test_diffwave_p3dfft ]; then
    alert "test_diffwave_p3dfft binary not found or not executable"
    exit 1
fi

if ! which mpiexec >/dev/null 2>/dev/null; then
    alert "Unable to find mpiexec; skipping test_diffwave_p3dfft"
    exit 77
fi

set -e # Fail on first error
mpiexec -np 1 ./test_diffwave_p3dfft
echo
mpiexec -np 2 ./test_diffwave_p3dfft
echo
mpiexec -np 4 ./test_diffwave_p3dfft
