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
banner()     { echo $banner_prefix: $*  ; }
run()        { echo $* ; $*             ; }
run_silent() { echo $* ; $* > /dev/null ; }

banner "Creating initial field to use for all tests"
run_silent ./channel_init --mms=0 --clobber "$tmpdir/initial.h5" \
                          --Nx=8 --Ny=16 --k=6 --htdelta=1 --Nz=6 \
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


banner "Equivalence of a field both with and without a restart"

# 1. Create an MMS field using channel_init at some Nx, Ny, k, Nz
# 2. Use channel_explicit --advance_nt=1 twice to advance to time steps passing through a restart file on disk.
# 3. Use channel_explicit --advance_nt=2 to advance two time steps.
# 4. Use h5diff to ensure /rho{,u,v,w,e} contents are identical between steps 2 and 3.


banner "Upsample/downsample both homogeneous directions"

# 1. Create an MMS field using channel_init at some Nx, Ny, k, Nz
# 2. Use channel_explicit --advance_nt=0 to upsample the field to 2*Nx, 3*Nz
# 3. Use channel_explicit --advance_nt=0 to downsample the field to Nx, Nz
# 4. Use h5diff to ensure /rho{,u,v,w,e} contents are identical between steps 1 and 3

banner "Upsample/downsample inhomogeneous direction order"

# 1. Create an MMS field using channel_init at some Nx, Ny, k, Nz
# 2. Use channel_explicit --advance_nt=0 to upsample the field to 2*k
# 3. Use channel_explicit --advance_nt=0 to downsample the field to k
# 4. Use h5diff to ensure /rho{,u,v,w,e} contents are equivalent between steps 1 and 3 to some tolerance


banner "Upsample/downsample inhomogeneous direction NDOF and htdelta"

# 1. Create an MMS field using channel_init at some Nx, Ny, k, htdelta, Nz
# 2. Use channel_explicit --advance_nt=0 to upsample the field to 2*Ny and 2*htdelta
# 3. Use channel_explicit --advance_nt=0 to downsample the field to Ny and htdelta
# 4. Use h5diff to ensure /rho{,u,v,w,e} contents are equivalent between steps 1 and 3 to some tolerance
