#!/bin/bash
# Produce a tab-separated table of filepath vs sorted times
# from HDF5 files containing time as a scalar in "/t"
set -eu

# Used to accumulate intermediate output
tmp1="$(mktemp)"
tmp2="$(mktemp)"
trap "rm -f \"$tmp1\" \"$tmp2\"" EXIT

# Extract time information from each file.
# h5dump could not make this much more awkward, could it?
# Especially since this does some initial newline hinky-ness.
while [[ $# -gt 0 ]]; do
    if [[ -f "$1" ]]; then
        h5dump -m %24.16f -o >(cat>>"$tmp1") -yd "/t" "$1" 1>/dev/null
        echo "$1" >> "$tmp2"
    else
        echo "$0: bad filepath \'$1\'" 1>&2
    fi
    shift
done

# Sort the combined results to standard output,
# computing time differences between adjacent rows
# to facilitate looking for "errant" samples
paste <(sed -e '/./,$!d' -e '$a\' "$tmp1") "$tmp2" \
    | sort -t$'\t' -n                              \
    | awk -F$'\t' '{$3=$1-t;print;t=$1}'           \
    | column -t
