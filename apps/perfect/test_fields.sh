#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_setup.sh"

# Ensure our restart-loading can handle sample laminar fields
banner "Restarting from current laminar restart file"
(
    : ${FIELDSDIR:=.}
    cd $testdir
    run ../channel "${FIELDSDIR}/laminar_k08.h5" --advance_nt=0
)

# Ensure our restart-loading routines remain backwards-compatible
banner "Restarting from legacy laminar restart file (r22804)"
(
    : ${FIELDSDIR:=.}
    cd $testdir
    run ../channel "${FIELDSDIR}/legacy_r22804.h5" --advance_nt=0
)
