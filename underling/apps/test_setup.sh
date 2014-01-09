#!/bin/bash
# Common functionality sourced from other test scripts

# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Check prerequisites and either warn and pass or die loudly
prereq_status=
for tool in egrep mpiexec
do
    if ! which $tool >/dev/null ; then
        echo "WARNING: Unable to find utility $tool" 1>&2
        prereq_status=0
    fi
done
for binary in ./underling_bench
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
banner() {
    echo -en "\n$banner_prefix:$1${METACASE:+ }"
    shift
    echo ${METACASE:+ (}${METACASE:-}${METACASE:+)}: "$@"
}
run()    { echo mpiexec -np 1        "$@" 1>&2 ; mpiexec -np 1        "$@"           ; }
prun()   { echo mpiexec -np ${NP:-1} "$@" 1>&2 ; mpiexec -np ${NP:-1} "$@"           ; }

# Create directory for scratch use
test -z "${TMPDIR-}" && export TMPDIR=.
testdir=`mktemp -d`

# Install teardown() function at exit unless TEST_UNDERLING_DEBUG is non-empty
declare -ir starttime=`date +%s`
teardown() {
    METACASE=
    banner "Tearing down"
    rm -rvf "$testdir"                     # Remove scratch directory
    test -z "`jobs -p`" || kill `jobs -p`  # Kill any lingering jobs

    declare -ir endtime=`date +%s`
    echo "execution took roughly $(expr $endtime - $starttime) seconds"
}
test -z "${TEST_UNDERLING_DEBUG-}" && trap "teardown" EXIT
