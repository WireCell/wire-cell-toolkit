#!/usr/bin/env bats
bats_load_library wct-bats.sh

@test "config export dump schema" {
    cd_tmp file

    mkdir -p raw log cfg uniq schema

    local cleanup="$(relative_path set-diff.jsonnet)"

    for libpath in $(find "$(blddir)" -name 'libWireCell*')
    do
        local libname="$(basename $libpath .so)"
        local name=${libname#"lib"}
        if [ "$name" == "WireCellApps" ] ; then
            extra=""
        else
            extra="-p $name"
        fi
        local out="raw/${name}.json"
        echo "[{type:\"ConfigSchema\",data:{filename: \"$out\"}}]" > \
             "cfg/${name}.jsonnet"

        wire-cell -l "log/${name}.log" -L debug -a ConfigSchema \
                  -p WireCellApps $extra "cfg/${name}.jsonnet"

    done

    # subtract the components provided by WireCellApps
    for one in raw/*.json
    do
        fname="$(basename $one)"
        if [ "$fname" == "WireCellApps.json" ] ; then
            continue;
        fi
        jsonnet -o "uniq/${fname}" \
                --tla-code-file small=raw/WireCellApps.json \
                --tla-code-file big=raw/${fname} \
                "$cleanup"
    done
    cp raw/WireCellApps.json uniq/

    # patch in missing schema 
    local patchup="$(relative_path add-missing-schema.jsonnet)"
    for one in uniq/*.json
    do
        pkg="$(basename $one .json)"
        jsonnet -o "schema/${pkg}.json" \
                --tla-code-file base="$one" -A pkg=$pkg $patchup
    done
}
