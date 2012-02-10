#!/bin/bash
# Common functionality sourced from multiple test scripts

# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Check prerequisites and either warn and pass or die loudly
prereq_status=
for tool in egrep cut mpiexec h5diff h5dump h5ls
do
    if ! which $tool >/dev/null ; then
        echo "WARNING: Unable to find utility $tool" 1>&2
        prereq_status=0
    fi
done
for binary in ./channel_init ./channel
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
banner() { echo; echo $banner_prefix${METACASE:+ (}${METACASE:-}${METACASE:+)}: "$@" ; }
run()    { echo mpiexec -np 1        "$@" ; mpiexec -np 1        "$@"                ; }
runq()   { echo mpiexec -np 1        "$@" ; mpiexec -np 1        "$@" > /dev/null    ; }
prun()   { echo mpiexec -np ${NP:-1} "$@" ; mpiexec -np ${NP:-1} "$@"                ; }
prunq()  { echo mpiexec -np ${NP:-1} "$@" ; mpiexec -np ${NP:-1} "$@" > /dev/null    ; }
differ() { echo h5diff "$@" ; h5diff "$@" || h5diff -rv "$@" ;}
differ_exclude() {
    h5diff_version_string=$(h5diff --version | tr -d '\n' | sed -e 's/^.*ersion  *//')
    h5diff_version_number=$(echo $h5diff_version_string | sed -e 's/\.//g')
    if test "$h5diff_version_number" -lt 186; then
        echo "WARN: Skipping 'h5diff $@' as $h5diff_version_string lacks --exclude-path"
    else
        echo h5diff "$@" ; h5diff "$@" || h5diff -rv "$@" ;
    fi
}

# Create directory for scratch use
test -z "${TMPDIR-}" && export TMPDIR=.
testdir=`mktemp -d`

# Install teardown() function at exit unless TEST_CHANNEL_DEBUG is non-empty
declare -ir starttime=`date +%s`
teardown() {
    METACASE=
    banner "Tearing down"
    rm -rvf "$testdir"                     # Remove scratch directory
    test -z "`jobs -p`" || kill `jobs -p`  # Kill any lingering jobs

    declare -ir endtime=`date +%s`
    echo "execution took roughly $(expr $endtime - $starttime) seconds"
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

banner "Building --exclude-paths for filtering samples"
datasets_bar=$(h5ls -f "$testdir/mms0.h5" | egrep '^/bar_' | cut "-d " -f 1)
exclude_datasets_bar=""
for dset in $datasets_bar
do
    exclude_datasets_bar="$exclude_datasets_bar --exclude-path=$dset"
done
exclude_datasets_bar=$(echo $exclude_datasets_bar | tr -d '\n' | tr -s ' ')
