#!/bin/bash
set -eu

# Run test template for explicit operators
OPER=--explicit
source "`dirname $0`/test_physical_template.sh"
