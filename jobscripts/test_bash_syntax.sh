#!/bin/bash
set -e

# Ensure 'file(1)' reports what we expect for a bash script
file "$BASH_SOURCE" | grep "Bourne-Again shell script" >/dev/null

# Check that all bash scripts have valid syntax
status=0
for entry in "*"; do
    if file "$entry" | grep "Bourne-Again shell script" >/dev/null; then
           (bash -n "$entry" && echo "$entry: syntax ok"   ) \
        || (status=1         && echo "$entry: syntax ERROR")
    fi
done
exit $status
