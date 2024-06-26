#!/usr/bin/bash
me="$(realpath "$BASH_SOURCE")"
mydir="$(dirname "$me")"

cd "$mydir"

for det in base pdsp uboone
do
    wirecell-util wires-volumes $det \
        | jsonnetfmt -n 4 --max-blank-lines 1 - \
                     > $mydir/mids/$det/variants/geometry.jsonnet
done
