#!/bin/bash
# Common functionality sourced from multiple test scripts

# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Check prerequisites and either warn and pass or die loudly
prereq_status=
for tool in mpiexec h5diff h5dump
do
    if ! which $tool >/dev/null ; then
        echo "WARNING: Unable to find utility $tool" 1>&2
        prereq_status=0
    fi
done
for binary in ./channel_init ./channel_explicit
do
    if [ ! -x $binary ]; then
        echo "ERROR: $binary not found or not executable" 1>&2
        prereq_status=1
    fi
done
if test x$prereq_status != x ; then
    echo `basename $0` ": unable to continue.  Exiting with status $prereq_status" 1>&2
    exit $prereq_status
fi

# Minimalistic command execution infrastructure
banner_prefix=`basename $0`
banner() { echo; echo $banner_prefix${METACASE+ (}${METACASE-}${METACASE+)}: "$@" ; }
run()    { echo mpiexec -np 1        "$@" ; mpiexec -np 1        "$@"             ; }
runq()   { echo mpiexec -np 1        "$@" ; mpiexec -np 1        "$@" > /dev/null ; }
prun()   { echo mpiexec -np ${NP:-1} "$@" ; mpiexec -np ${NP:-1} "$@"             ; }
prunq()  { echo mpiexec -np ${NP:-1} "$@" ; mpiexec -np ${NP:-1} "$@" > /dev/null ; }
differ() { echo h5diff "$@" ; h5diff "$@" || h5diff -rv "$@" ;}

# Create directory for scratch use
testdir=`mktemp -d`

yes >/dev/null &
# Install teardown() function at exit unless TEST_CHANNEL_DEBUG is non-empty
teardown() {
    METACASE=
    banner "Tearing down"
    rm -rvf "$testdir"                     # Remove scratch directory
    test -z "`jobs -p`" || kill `jobs -p`  # Kill any lingering jobs
}
test -z "${TEST_CHANNEL_DEBUG-}" && trap "teardown" EXIT

banner "Creating initial field to use for tests"
declare -ir Nx=4
declare -ir Ny=12
declare -ir k=6
declare -ir htdelta=1
declare -ir Nz=6
runq ./channel_init "$testdir/mms0.h5" --mms=0                           \
                    --Nx=$Nx --Ny=$Ny --k=$k --htdelta=$htdelta --Nz=$Nz
chmod +r "$testdir/mms0.h5"
