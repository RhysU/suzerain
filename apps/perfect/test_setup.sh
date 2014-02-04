#!/bin/bash
# Common functionality sourced from multiple test scripts

# Initialize test infrastructure
source "`dirname $0`/../test_infrastructure.sh"

# Check prerequisites and die loudly if the tools we need aren't available
prereq_status=
for binary in ./perfect_initial ./perfect_advance
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


# $(testdir) has already been created by test_infrastructure.sh for scratch use

# We save the initial wisdom in $testdir for possible reuse (ticket #2515)
banner "Creating initial field to use for tests"
declare -r  Re=100
declare -r  Ma=1.15
declare -ir Nx=4
declare -ir Ny=12
declare -ir k=6
declare -ir htdelta=1
declare -ir Nz=6
runq ./perfect_initial "$testdir/mms0.h5" --mms=0 --Re=$Re --Ma=$Ma         \
                       --Nx=$Nx --Ny=$Ny --k=$k --htdelta=$htdelta --Nz=$Nz \
                       ${DECOMP:-} "--plan_wisdom=$testdir/wisdom.init"
chmod +r "$testdir/mms0.h5"


banner "Checking zero-zero and Nyquist wavenumbers are strictly real-valued"
(
    cd $testdir

    # Build a single h5dump invocation dumping the relevant wavenumbers
    declare -a realmodes=("0,0,0")
    declare -a cmd=("${H5DUMP}")
    [ "$((Nx%2))" -eq 0 ]                && realmodes+=("0,$((Nx/2)),0")
    [ "$((Nz%2))" -eq 0 ]                && realmodes+=("$((Nz/2)),0,0")
    [ "$((Nx%2))" -eq 0 -a "$((Nz%2))" ] && realmodes+=("$((Nz/2)),$((Nx/2)),0")
    for field in rho rho_u rho_v rho_w rho_E
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

    # Ensure that the entire fourth column (imag component) is identically zero
    awk $awk_h5dump_dataonly realmodes | awk '$4 != 0.0 { exit 1 }' \
        || (cat realmodes && false)
)


banner "Building --exclude-paths for filtering samples"
datasets_bar=$(${H5LS} -f "$testdir/mms0.h5" | ${GREP} '^/bar_' | cut "-d " -f 1)
exclude_datasets_bar=""
for dset in $datasets_bar
do
    exclude_datasets_bar="$exclude_datasets_bar --exclude-path=$dset"
done
exclude_datasets_bar=$(echo $exclude_datasets_bar | tr -d '\n' | tr -s ' ')
