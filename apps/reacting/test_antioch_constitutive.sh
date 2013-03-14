#!/bin/bash
set -eu

# Initialize test infrastructure
source "`dirname $0`/../test_infrastructure.sh"

# Run all antioch tests
banner "Test antioch_constitutive class capabilities"
(
    : ${ANTIOCH_DATA_DIR:=.}
    cd $testdir # created by test_infrastructure.sh
    ../test_antioch_constitutive "${ANTIOCH_DATA_DIR}/air_5sp.xml"
)
