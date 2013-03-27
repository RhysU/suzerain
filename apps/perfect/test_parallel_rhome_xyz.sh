#!/bin/bash
set -eu

# Run test template for implicit operators
OPER=--implicit=rhome_xyz
source "`dirname $0`/test_parallel_template.sh"
