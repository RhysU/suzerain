#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Fail on first error or undefined value
set -eu

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
differ() { echo h5diff "$@" ; h5diff "$@" || h5diff -r "$@" ;}

banner "Creating initial field to use for all tests"
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

banner "Idempotence of restarting without time advancement"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0
    differ initial.h5 a0.h5
)

banner "Equivalence of a field both with and without a restart"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=1
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=1
    runq ../channel_explicit initial.h5 --desttemplate "c#.h5" --advance_nt=2
    differ --use-system-epsilon b0.h5 c0.h5
)

banner "Upsample/downsample both homogeneous directions"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                        --Nx=$((2*$Nx)) --Nz=$((3*$Nz))
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=0 \
                                        --Nx=$((  $Nx)) --Nz=$((  $Nz))
    differ initial.h5 b0.h5
)

banner "Upsample/downsample inhomogeneous direction order"
(
    cd $tmpdir
    runq ../channel_explicit initial.h5 --desttemplate "a#.h5" --advance_nt=0 \
                                        --k=$(($k+1))
    runq ../channel_explicit a0.h5      --desttemplate "b#.h5" --advance_nt=0 \
                                        --k=$(($k  ))
    # Chosen tolerances are wholly empirical and represent nothing deep
    differ --delta=5e-5 initial.h5 b0.h5 /rho
    differ --delta=3e-4 initial.h5 b0.h5 /rhou
    differ --delta=5e-5 initial.h5 b0.h5 /rhov
    differ --delta=6e-5 initial.h5 b0.h5 /rhow
    differ --delta=7e-4 initial.h5 b0.h5 /rhoe
)

banner "Upsample/downsample inhomogeneous direction NDOF and htdelta"

# 1. Create an MMS field using channel_init at some Nx, Ny, k, htdelta, Nz
# 2. Use channel_explicit --advance_nt=0 to upsample the field to 2*Ny and 2*htdelta
# 3. Use channel_explicit --advance_nt=0 to downsample the field to Ny and htdelta
# 4. Use h5diff to ensure /rho{,u,v,w,e} contents are equivalent between steps 1 and 3 to some tolerance
