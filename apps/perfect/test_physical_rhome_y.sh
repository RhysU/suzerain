#!/bin/bash
set -eu

# Run test template for implicit operators in Y direction
OPER=--implicit=rhome_y
source "`dirname $0`/test_physical_template.sh"
