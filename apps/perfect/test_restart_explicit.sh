#!/bin/bash
set -eu

# Run test template for explicit operators
OPER=--explicit
source "`dirname $0`/test_restart_template.sh"
