#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# Ensure our restart-loading can handle sample laminar fields
banner "Restarting from current laminar restart file"
(
    : ${FIELDSDIR:=.}
    cd $testdir
    run ../perfect_advance ${DECOMP:-} "${FIELDSDIR}/channel_k08.h5" \
        --advance_nt=0
)

# Ensure our restart-loading routines remain backwards-compatible
for name in legacy_r{22804,38808}.h5;
do
    banner "Restarting from legacy laminar restart file ($name)"
    (
        : ${FIELDSDIR:=.}
        cd $testdir
        run ../perfect_advance ${DECOMP:-} "${FIELDSDIR}/$name" \
            --advance_nt=0
    )
done
