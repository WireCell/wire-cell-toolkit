#!/usr/bin/env bash

# in-place generator for mids.jsonnet to save some typing.

me="$(realpath "$BASH_SOURCE")"
mydir="$(dirname "$me")"

cd "$mydir"

out="mids.jsonnet"

cat <<EOF > "$out"
// This file is generated, do not edit.
// Generated on $(date) by $USER on $(hostname --fqdn)
// with $me
{
EOF

# this loop should keep path relative
for one in mids/*/mids.jsonnet
do
    det=$(basename $(dirname $one))
    echo "    $det : import \"$one\"," >> "$out"
done
echo "}" >> "$out"

    
