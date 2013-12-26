#!/bin/bash
# Use latexmk to build one writeup without the full Suzerain build system.
# For example, 'writeup baseflow.tex' should build 'baseflow.pdf'.
# Be aware that will tend to be exceedingly brittle.
set -eu
target=$1
shift
SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
TEXINPUTS="${TEXINPUTS:+$TEXINPUTS:}${SCRIPTDIR}:${SCRIPTDIR}/.." \
    latexmk                                                       \
    "-auxdir=${TMPDIR-/tmp}"                                      \
    -dvi-                                                         \
    -ps-                                                          \
    -pdf                                                          \
    -recorder                                                     \
    "$@"                                                          \
    -pdflatex="pdflatex --shell-escape %O %S" "$target"
