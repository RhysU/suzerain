#!/bin/bash
set -eu

# Run test template for implicit operators
OPER=--implicit
source "`dirname $0`/test_restart_template.sh"
