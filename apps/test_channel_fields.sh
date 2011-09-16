#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/test_channel_setup.sh"

# Ensure our restart-loading routines remain backwards-compatible
banner "Restarting from legacy laminar restart file"
(
    : ${FIELDSDIR:=.}
    cd $testdir
    runq ../channel_explicit "${FIELDSDIR}/laminar_k08".h5 --advance_nt=0
)
