#!/bin/bash
# Common functionality sourced from multiple test scripts

# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

# Check prerequisites and die loudly if the tools we need aren't available
prereq_status=
for tool in column cut mktemp mpiexec tail
do
    if ! which $tool >/dev/null 2>/dev/null; then
        echo "ERROR: Unable to find utility $tool" 1>&2
        prereq_status=1
    fi
done
if test x$prereq_status != x ; then
    echo `basename $0` ": unable to continue.  Exiting with status $prereq_status" 1>&2
    exit $prereq_status
fi

# Create directory for scratch use
test -z "${TMPDIR-}" && export TMPDIR=.
testdir="$(readlink -f "$(mktemp -d "$TMPDIR/testdir.XXXXXX")")"

# Install teardown() function at exit unless TEST_DEBUG is non-empty
declare -ir starttime=`date +%s`
teardown() {
    METACASE=
    banner "Tearing down"
    rm -rvf "$testdir"                     # Remove scratch directory
    test -z "`jobs -p`" || kill `jobs -p`  # Kill any lingering jobs

    declare -ir endtime=`date +%s`
    echo "execution took roughly $(expr $endtime - $starttime) seconds"
}
test -z "${TEST_DEBUG-}" && trap "teardown" EXIT

# Minimalistic command execution infrastructure

alert() {
    # Requires AM_TESTS_FD_REDIRECT in Makefile.am
    if test -t 9; then
        echo WARN: "$@" >&9
    fi
    echo WARN: "$@"
}

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
    local h5diff_version_string=$(${H5DIFF} --version | tr -d '\n' | sed -e 's/^.*ersion  *//' -e 's/-patch.*$//')
    local h5diff_version_number=$(echo $h5diff_version_string | sed -e 's/\.//g')
    if test "$h5diff_version_number" -lt 186; then
        alert "Skipping portions of test as ${H5DIFF} $h5diff_version_string lacks required --exclude-path"
        exit 77 # See http://www.gnu.org/software/automake/manual/html_node/Scripts_002dbased-Testsuites.html
    fi
    local outfile=$(mktemp "$testdir/differ.XXXXXX")
    local prefix="--exclude-path /metadata_generated"
    echo ${H5DIFF} --verbose --report $prefix "$@" 2>&1 | tee -a $outfile
    echo                                           2>&1 |      >>$outfile
    # tail not cat because h5diff echos its invocation arguments
    # embedded awk script used to add a ratio column and pretty up the output
    ${H5DIFF} --verbose --report $prefix "$@"      2>&1        >>$outfile || (tail -n +2 $outfile | ${GREP} -v "^0 differences found$" | ${AWK} -f <(cat - <<-'HERE'
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
                if (b[2]+0 == 0) {
                    # Output nothing otherwise division by zero occurs
                    # Outputting "NaN" is unsatisfactory as Inf might result
                } else {
                    printf OFMT, b[1]/b[2] | aligner
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
