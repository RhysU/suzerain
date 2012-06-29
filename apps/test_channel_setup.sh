#!/bin/bash
# Common functionality sourced from multiple test scripts

# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Check prerequisites and die loudly if the tools we need aren't available
prereq_status=
for tool in egrep cut mpiexec h5diff h5dump h5ls
do
    if ! which $tool >/dev/null ; then
        echo "ERROR: Unable to find utility $tool" 1>&2
        prereq_status=1
    fi
done
for binary in ./channel_init ./channel
do
    if [ ! -x $binary ]; then
        echo "ERROR: $binary not found or not executable" 1>&2
        prereq_status=2
    fi
done
if test x$prereq_status != x ; then
    echo `basename $0` ": unable to continue.  Exiting with status $prereq_status" 1>&2
    exit $prereq_status
fi

# Create directory for scratch use
test -z "${TMPDIR-}" && export TMPDIR=.
testdir="$(readlink -f "$(mktemp -d)")"

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

# Minimalistic command execution infrastructure
banner_prefix=`basename $0`
banner() {
    echo
    msg="$banner_prefix${METACASE:+ (}${METACASE:-}${METACASE:+)}: $@"
    printf "%$(expr length "$msg")s\n"|tr ' ' '#'
    echo $msg
    printf "%$(expr length "$msg")s\n"|tr ' ' '#'
    echo
}
run()    { echo mpiexec -np 1        "$@" ; mpiexec -np 1        "$@"             ; echo; }
runq()   { echo mpiexec -np 1        "$@" ; mpiexec -np 1        "$@" > /dev/null ; echo; }
prun()   { echo mpiexec -np ${NP:-1} "$@" ; mpiexec -np ${NP:-1} "$@"             ; echo; }
prunq()  { echo mpiexec -np ${NP:-1} "$@" ; mpiexec -np ${NP:-1} "$@" > /dev/null ; echo; }
differ() {
    outfile=`mktemp --tmpdir="$testdir"`
    echo h5diff "$@"  2>&1 | tee -a $outfile
    echo              2>&1 |      >>$outfile
    # tail not cat because h5diff echos its invocation arguments
    # embedded awk script used to add a ratio column and pretty up the output
    h5diff -r "$@"    2>&1        >>$outfile || (tail -n +2 $outfile | egrep -v "^0 differences found$" | awk -f <(cat - <<-'HERE'
            BEGIN { OFMT=" %+14.8g"; aligner="column -t" }
            {
                sub("[[:space:]]*$", "")
            }
            $0 ~ /[[:space:]]+difference([[:space:]]+relative)?$/ {
                print $0"    ratio"
                getline
                print $0"---------"
                augment = 1
                next
            }
            $0 ~ /^[[:space:]]*[[:digit:]]+[[:space:]]+differences?[[:space:]]+found$/ {
                close(aligner)
                print
                augment = 0
                next;
            }
            augment == 1 {
                split($0,   a, /\][[:space:]]+/)
                printf "%s ]", a[1] | aligner
                n = split(a[2], b)
                for (i = 1; i <= n; ++i) {
                    printf OFMT, b[i] | aligner
                }
                if (b[2] != 0) {
                    printf OFMT, b[1]/b[2] | aligner
                } else {
                    printf " NaN" | aligner
                }
                printf "\n"| aligner
                next
            }
            {
                print
            }
HERE
    ) && false)
}
differ_exclude() {
    h5diff_version_string=$(h5diff --version | tr -d '\n' | sed -e 's/^.*ersion  *//')
    h5diff_version_number=$(echo $h5diff_version_string | sed -e 's/\.//g')
    if test "$h5diff_version_number" -lt 186; then
        echo "WARN: Skipping 'h5diff $@' as $h5diff_version_string lacks --exclude-path\n"
    else
        differ "$@"
    fi
}


banner "Creating initial field to use for tests"
declare -r  Re=100
declare -r  Ma=1.15
declare -ir Nx=4
declare -ir Ny=12
declare -ir k=6
declare -ir htdelta=1
declare -ir Nz=6
runq ./channel_init "$testdir/mms0.h5" --mms=0 --Re=$Re --Ma=$Ma         \
                    --Nx=$Nx --Ny=$Ny --k=$k --htdelta=$htdelta --Nz=$Nz
chmod +r "$testdir/mms0.h5"


banner "Checking zero-zero and Nyquist wavenumbers are strictly real-valued"
(
    cd $testdir

    # Build a single h5dump invocation dumping the relevant wavenumbers
    declare -a realmodes=("0,0,0")
    declare -a cmd=('h5dump')
    [ "$((Nx%2))" -eq 0 ]                && realmodes+=("0,$((Nx/2)),0")
    [ "$((Nz%2))" -eq 0 ]                && realmodes+=("$((Nz/2)),0,0")
    [ "$((Nx%2))" -eq 0 -a "$((Nz%2))" ] && realmodes+=("$((Nz/2)),$((Nx/2)),0")
    for field in rho rhou rhov rhow rhoe
    do
        for realmode in "${realmodes[@]}"
        do
            cmd+=(-d "/$field[$realmode;;1,1,$Ny]")
        done
    done
    cmd+=(-w 8 -m %22.14g mms0.h5)
    echo ${cmd[*]}
    ${cmd[*]} > realmodes 2>&1 \
        || (cat realmodes && false)

    # awk program keeping only the h5dump lines between 'DATA {' and '}'
    declare -r awk_h5dump_dataonly='                                        \
        $0 ~ /^[[:space:]]*DATA[[:space:]]*{[[:space:]]*$/ {keep = 1; next} \
        $0 ~ /^[[:space:]]*}[[:space:]]*$/                 {keep = 0}       \
        keep == 1 { print $0 }'

    # Ensure that the entire fourth column (the imaginary component) is identically zero
    awk $awk_h5dump_dataonly realmodes | awk '$4 != 0.0 { exit 1 }' \
        || (cat realmodes && false)
)


banner "Building --exclude-paths for filtering samples"
datasets_bar=$(h5ls -f "$testdir/mms0.h5" | egrep '^/bar_' | cut "-d " -f 1)
exclude_datasets_bar=""
for dset in $datasets_bar
do
    exclude_datasets_bar="$exclude_datasets_bar --exclude-path=$dset"
done
exclude_datasets_bar=$(echo $exclude_datasets_bar | tr -d '\n' | tr -s ' ')
