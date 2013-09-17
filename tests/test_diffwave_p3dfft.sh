#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x test_diffwave_p3dfft ]; then
    echo "test_diffwave_p3dfft binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping test_diffwave_p3dfft"
    exit 0
fi

set -e # Fail on first error
mpiexec -np 1 ./test_diffwave_p3dfft
echo
mpiexec -np 2 ./test_diffwave_p3dfft
echo
mpiexec -np 4 ./test_diffwave_p3dfft
