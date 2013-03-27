#!/bin/bash
set -eu

# Run test template for implicit operators in X, Y, and Z directions
OPER=--implicit=rhome_xyz
source "`dirname $0`/test_parallel_template.sh"
