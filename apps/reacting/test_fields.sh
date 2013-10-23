#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# TODO: Switch back to implicit once that capability is resurrected

# Ensure our restart-loading can handle sample laminar fields
banner "Restarting from current laminar restart file"
(
    : ${FIELDSDIR:=.}
    cd $testdir
    run ../reacting_advance --explicit "${FIELDSDIR}/channel_k08.h5" --advance_nt=0
)

# Ensure our restart-loading routines remain backwards-compatible
banner "Restarting from legacy laminar restart file (r22804)"
(
    : ${FIELDSDIR:=.}
    cd $testdir
    run ../reacting_advance --explicit "${FIELDSDIR}/legacy_r22804.h5" --advance_nt=0
)
