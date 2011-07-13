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

# Process command line options
# TODO Allow -np to specify a number of MPI ranks to use
# TODO Accept -v as a verbosity flag
# TODO Pass all other arguments through verbatim

# Create temporary directory and clean it up on exit
tmpdir=`mktemp -d`
trap "rm -rf $tmpdir" EXIT

# Minimalistic command execution infrastructure
banner_prefix=`basename $0`
banner() { echo; echo $banner_prefix: "$@" ; }
run()    { echo "$@" ; "$@"                ; }
runq()   { echo "$@" ; "$@" > /dev/null    ; }
differ() { echo h5diff "$@" ; h5diff "$@" || h5diff -rv "$@" ;}

banner "Creating initial field to use for tests"
runq ./channel_init "$tmpdir/initial.h5"                            \
                    --mms=0 --Nx=4 --Ny=12 --k=6 --htdelta=1 --Nz=6 \
                    $* # Incoming script arguments override

# Slurp grid details from the restart into integer variables
# This account for any overrides present on the channel_init line just above
function read_restart() {
    h5dump -y -d $1 -o "$tmpdir/read_restart" "$tmpdir/initial.h5" >/dev/null
    cat "$tmpdir/read_restart"
}
declare -ir Nx=$(read_restart Nx)
declare -ir Ny=$(read_restart Ny)
declare -ir k=$(read_restart k)
declare -ir htdelta=$(read_restart htdelta)
declare -ir Nz=$(read_restart Nz)
