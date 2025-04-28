#!/usr/bin/env bats

# This is derived from various parts of:

base_url="https://github.com/HaiwangYu/wcp-porting-img"
raw_url="$base_url/raw/refs/heads/main"

bats_load_library wct-bats.sh

# Process a log from one of the tests into a smaller "digest" log and check if
# it changed if a historical version exists.
do_log_digest () {
    local log="$1"; shift              # the log file name
    debug "Checking log: $log"
    test -s "$log"

    local ign="${1:-downloads}"
    local dig="${log%.log}.dig"

    sort $log | grep ' D ' | grep '\[' | egrep -v "$ign" | while read line
    do
        printf "%s\n" "${line:15}"
    done > $dig

    debug "Checking log digest: $dig"
    test -s $dig

    old_dig=$(historical_files --allow-missing --last $(current-test-name)/$dig)    
    if [ -n "$old_dig" ] ; then
        diff $dig $old_dig
    fi
    saveout -c history $dig
}    


do_prep () {
    name="$1"
    cd_tmp file
    mkdir -p $name
    cd $name
}

@test "porting qlport" {
    local name=qlport
    do_prep $name

    local url="$raw_url/$name/rootfiles/nuselEval_5384_130_6501.root"

    local dat="$(download_file "$url")"
    local cfg="$(relative_path test-porting/$name/main.jsonnet)"
    local bee="$name.zip"
    local log="$name.log"
    local dig="$name.dig"
    local dag="$name.pdf"

    test -n "$dat"

    run_idempotently -s "$cfg" -t "$dag" -- \
                     wirecell-pgraph dotify $cfg $dag \
                     -A kind=both -A "infiles=$dat" -A "beezip=$bee"


    run_idempotently -s "$cfg" -s "$dat" -t "$bee" -t "$log" -- \
                     wire-cell -l "$log" -L debug \
                     -A kind=both -A "infiles=$dat" -A "beezip=$bee" "$cfg" 
    do_log_digest $log $dat

    for zip in qlport.zip
    do
        file_larger_than $zip 22
    done
}


@test "porting pdhd" {

    local name=pdhd
    do_prep $name

    local gitref=b5f7f21e1ca853e29d746ae9044ac79c885956b0
    local indir=$(download_git_subdir --ref $gitref "$base_url" pdhd/1event)
    debug "PDHD input from $indir"
    test -n "$indir"
    test -d "$indir"
    test -s "$indir/clusters-apa-apa0.tar.gz"

    local cfg="$(relative_path test-porting/$name/main.jsonnet)"
    local bee="mabc-all-apa.zip"
    local log="$name.log"
    local dig="$name.dig"

    run_idempotently -s "$cfg" -s "$dat" -t "$bee" -t "$log" -- \
                     wire-cell -l "$log" -L debug \
                     -A input=$indir "$cfg" 

    # Note, "empty" zips have finite size, but none should 0 bytes.
    for zip in *.zip ; do
        test -s $zip
    done

    for zip in mabc-all-apa.zip mabc-apa0-face0.zip mabc-apa1-face1.zip mabc-apa2-face0.zip mabc-apa3-face1.zip
    do
        file_larger_than $zip 22
    done

    do_log_digest $log $indir
}

@test "porting fgval" {

    local name=fgval
    do_prep $name

    local run=5384
    local sub=130
    local evt=6501

    local url="$raw_url/$name/result_${run}_${sub}_${evt}.root"
    local cfg1="$(relative_path test-porting/$name/stage1.jsonnet)"

    local dat="$(download_file "$url")"
    test -s "$dat"
    debug "Input file: $dat"

    for what in live dead ; do
        mkdir -p $what

        local out="$what/clusters.npz"
        local log="$what/clusters.log"

        run_idempotently -s "$cfg1" -s "$dat" -t "$out" -t "$log" -- \
                     wire-cell -l "$log" -L debug -A iname="$dat" -A oname="$out" -A kind="$what" "$cfg1"
        do_log_digest "$log"
    done

    local cfg2="$(relative_path test-porting/$name/stage2.jsonnet)"
    local log="clustering.log"
    local bee="clustering.zip"
    local out="tensor-apa-uboone.tar.gz"

    run_idempotently -s "$cfg2" -s "live/clusters.npz" -s "dead/clusters.npz" -t "$bee" -t "$log" -t "$out" -- \
                     wire-cell -l "$log" -L debug \
                     -A "active_clusters=live/clusters.npz" \
                     -A "masked_clusters=dead/clusters.npz" \
                     -A "bee_zip=$bee" \
                     -A initial_index=0 \
                     -A initial_runNo=$run \
                     -A initial_subRunNo=$sub \
                     -A initial_eventNo=$evt \
                     "$cfg2"
    do_log_digest "$log" npz


}
