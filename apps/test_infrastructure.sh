#!/bin/bash
# Common functionality sourced from multiple test scripts
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
if test x$prereq_status != x ; then
    echo `basename $0` ": unable to continue.  Exiting with status $prereq_status" 1>&2
    exit $prereq_status
fi

# Create directory for scratch use
test -z "${TMPDIR-}" && export TMPDIR=.
testdir="$(readlink -f "$(mktemp -d)")"

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
