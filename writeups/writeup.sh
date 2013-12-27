#!/bin/bash
# Use latexmk to build one writeup without the full Suzerain build system.
# For example, 'writeup.sh baseflow.tex' should build 'baseflow.pdf'.
# Notably, 'writeup.sh baseflow.tex -pvc' starts a continuous previewer.
# Be aware that this hack may be exceedingly brittle.
set -eu
target=$1
shift
SCRIPTDIR="$( cd "$( echo "${BASH_SOURCE[0]%/*}" )"; pwd )"
TEXINPUTS=":${TEXINPUTS:+$TEXINPUTS:}${SCRIPTDIR}:${SCRIPTDIR}/.." \
BIBINPUTS=":${BIBINPUTS:+$BIBINPUTS:}${SCRIPTDIR}:${SCRIPTDIR}/.." \
    latexmk -dvi- -ps- -pdf -recorder                              \
    "$@"                                                           \
    -pdflatex="pdflatex --shell-escape %O %S"                      \
    "$SCRIPTDIR/$target"
