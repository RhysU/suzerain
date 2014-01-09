#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x test_underling ]; then
    echo "test_underling binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping test_underling"
    exit 0
fi

for np in 1 2 3 4
do
    cmd="mpiexec -np $np ./test_underling_fftw $@"
    echo $cmd
    $cmd
done
