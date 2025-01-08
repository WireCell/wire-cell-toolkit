#!/usr/bin/env bats

bats_load_library wct-bats.sh

@test "logging" {
    debug "debug"
#    log "log"
    info "info"
    warn "warn"
    error "error"
    fatal "fatal"
#    die "die"
}
    

@test "index of" {
    local a="$(index_of a a b c)"
    [[ "$a" = "0" ]]
    local b="$(index_of b a b c)"
    [[ "$b" = "1" ]]
    local abc=( a "b 1" "c d e f" )
    local c="$(index_of "c d e f" "${abc[@]}")"
    [[ "$c" = "2" ]]
}

@test "check check" {
    check true
    local failed=no
    check false || failed=yes
    [[ "$failed" = "yes" ]]
}

@test "using jq" {
    skip_if_missing jq

    cd_tmp
    output="$(echo '{"a":42}' | jq '.a')"
    [[ "$output" = "42" ]]
    [[ -n "$(jq -h | head -1 | grep 'jq - commandline JSON processor')" ]]
    echo '{"a":42}'  > junk.txt
    [[ "$(jq '.a' < junk.txt)" = "42" ]]
}

@test "check usepkg" {
    usepkg cfg
    [[ -n "$cfg_src" ]]
    [[ -d "$cfg_src" ]]

    cd_tmp
    usepkg util
    [[ -n "$util_src" ]]
    [[ -d "$util_src" ]]
}

@test "resolve file in src" {
    local got="$(resolve_file test/test/test_wct_bats.bats)"
    [[ -n "$got" ]]
    [[ -n "$BATS_TEST_FILENAME" ]]
    a="$(cat $BATS_TEST_FILENAME | md5sum | awk '{print $1}')"
    b="$(cat $got | md5sum | awk '{print $1}')"
    debug "md5: $a $BATS_TEST_FILENAME"
    debug "md5: $b $got"
    [[ "$a" = "$b" ]]
}

@test "tojson" {
    check tojson foo=42 bar="baz=quax quax"
    [[ -n "$(echo $output | grep '"foo": "42"')" ]]
    [[ -n "$(echo $output | grep '"bar": "baz=quax quax"')" ]]
}

@test "wcb env" {
    # disable this as it has become common to define PREFIX in direnv
    # [[ -z "$PREFIX" ]]
    # local prefix=$(wcb_env_value PREFIX)
    # [[ -z "$(echo $prefix | grep '"')" ]]
    # [[ -z "$PREFIX" ]]

    # only one var
    check wcb_env_vars VERSION
    [[ -n "$(echo $output | grep VERSION)" ]]
    [[ -z "$(echo $output | grep SUBDIRS)" ]]

    # all vars
    check wcb_env_vars 
    [[ -n "$(echo $output | grep VERSION)" ]]
    [[ -n "$(echo $output | grep SUBDIRS)" ]]

    # specific environment 
    [[ -z "$VERSION" ]]
    wcb_env VERSION
    [[ -n "$VERSION" ]]
    [[ -z "$SUBDIRS" ]]

    # everything
    wcb_env
    [[ -n "$VERSION" ]]
    [[ -n "$SUBDIRS" ]]

}

    # disable this as it has become common to define PREFIX in direnv
# @test "wcb env in temp dir" {
#     cd_tmp
#     [[ -z "$PREFIX" ]]
#     local prefix=$(wcb_env_value PREFIX)
#     [[ -z "$(echo $prefix | grep '"')" ]]
#     [[ -z "$PREFIX" ]]
# }


@test "have a test data file or skip" {
    local name="depos/many.tar.bz2"

    skip_if_no_input "$name"

    local got=$(input_file $name)
    [[ -s $got ]]
}

@test "relative path" {
    local got="$(relative_path check_batsio.bats)"
    echo "$got"
    [[ -n "$got" ]]
    [[ -f "$got" ]]
    [[ -n "$(echo $got | grep '^/')" ]]
}

@test "tmp dirs" {
    if [ -n "$WCT_BATS_TMPDIR" ] ; then
        [[ "$(tmpdir)" = "$WCT_BATS_TMPDIR" ]]
        return
    fi

    [[ -n "$BATS_RUN_TMPDIR" ]]
    [[ -n "$BATS_FILE_TMPDIR" ]]
    [[ -n "$BATS_TEST_TMPDIR" ]]
    [[ -n "$BATS_TEST_TMPDIR" ]]

    [[ -n "$(tmpdir run)" ]]
    [[ -n "$(tmpdir file)" ]]
    [[ -n "$(tmpdir test)" ]]
    [[ -n "$(tmpdir giraffea)" ]]

    [[ "$(tmpdir run)" = "$BATS_RUN_TMPDIR" ]]
    [[ "$(tmpdir file)" = "$BATS_FILE_TMPDIR" ]]
    [[ "$(tmpdir test)" = "$BATS_TEST_TMPDIR" ]]
    [[ "$(tmpdir giraffea)" = "$BATS_TEST_TMPDIR" ]]

    cd_tmp 
    [[ "$(pwd)" = "$BATS_TEST_TMPDIR" ]]
    cd_tmp file
    [[ "$(pwd)" = "$BATS_FILE_TMPDIR" ]]
}


@test "bats run command env var" {
    check env
    echo "run command: $BATS_RUN_COMMAND"
    [[ -n "$BATS_RUN_COMMAND" ]]
    # why is it not in env?
    [[ -z "$(env |grep BATS_RUN_COMMAND)" ]]
    export BATS_RUN_COMMAND
    [[ -n "$(env |grep BATS_RUN_COMMAND)" ]]
}
    
@test "topdir" {
    local t;
    t=$(topdir)
    run diff "${BASH_SOURCE[0]}" "$t/test/test/test_wct_bats.bats"
}

@test "pathlist resolution" {
    local t=$(topdir)
    declare -a got
    got=( $(resolve_pathlist $t/gen:$t/test $t/sigproc:$t/apps $t/doesnotexist) )
    # yell "got: ${#got[@]}: |${got[@]}|"
    [[ ${#got[@]} -eq 4 ]]
    [[ ${got[0]} = $t/gen ]]
    [[ ${got[1]} = $t/test ]]
    [[ ${got[2]} = $t/sigproc ]]
    [[ ${got[3]} = $t/apps ]]
}

@test "category resolution" {
    
    local category="testing"
    local test_data="$(blddir)/tests"
    local release=0.0.0
    local dirty="${release}-dirty"

    for one in $release $dirty
    do
        dir=$test_data/$category/$one
        rm -rf $dir
        mkdir -p $dir
        touch $dir/junk.txt
    done

    declare -a got=( $(find_category_version_paths -c testing) )
    [[ ${#got[@]} -eq 1 ]]
    [[ ${got[0]} = "$test_data/testing/$release" ]]

    got=( $(find_category_version_paths -c testing -d) )
    [[ ${#got[@]} -eq 2 ]]
    [[ ${got[0]} = "$test_data/testing/$release" ]]
    [[ ${got[1]} = "$test_data/testing/$dirty" ]]

    got=( $(find_category_version_paths -c testing -d) )
    [[ ${#got[@]} -eq 2 ]]
    [[ ${got[0]} = "$test_data/testing/$release" ]]
    [[ ${got[1]} = "$test_data/testing/$dirty" ]]
    

    declare -a got3=( $(find_paths junk.txt ${got[*]}) )
    [[ ${#got3[@]} -eq 2 ]]
    [[ "$(basename ${got3[0]})" = "junk.txt" ]]
    [[ "$(basename ${got3[1]})" = "junk.txt" ]]


    vers=( $(for n in ${got[*]}; do basename $n; done | sort -u | tail -1 ) )
    [[ "$vers" = "0.0.0-dirty" ]]
    

}

@test "run idempotently" {
    cd_tmp

    for f in file1.txt file2.txt
    do
        rm -f $f
        # unconditinally run
        run_idempotently -- touch $f
        [[ -f $f ]]
    done
    
    rm -f file3.txt
    run_idempotently -s file1.txt -t file3.txt -- cp file1.txt file3.txt
    [[ -f file3.txt ]]

    sync     # can run fast enough that file3.txt isn't fully created?
    touch timestamp.txt   # make this just to mark the time for later.

    # this should not update file3.txt
    run_idempotently -s file1.txt -t file3.txt -- cp file1.txt file3.txt
    [[ -f file3.txt ]]
    [[ file3.txt -ot timestamp.txt ]]
    [[ timestamp.txt -nt file3.txt ]]

    # this should update file3
    run_idempotently -s timestamp.txt -t file3.txt -- cp file1.txt file3.txt
    [[ -f file3.txt ]]
    [[ file1.txt -ot file3.txt ]]
    [[ file3.txt -nt file1.txt ]]

}

@test "historical versions" {
    # IFS=" " read -r -a versions <<< "$(historical_versions)"
    local verlines=$(historical_versions)
    local versions=( $verlines ) # split
    # yell "versions: ${versions[@]}"
    # yell "versions[0]: ${versions[0]}"
    [[ "${#versions[@]}" -gt 0 ]]
    [ -n "${versions[0]}" ]
}
@test "historical files" {
    local hf
    hf=( $(historical_files test-addnoise/frames-adc.tar.gz) )
    debug "hf all: ${#hf[@]}: ${hf[*]}"
    hf=( $(historical_files --last 2 test-addnoise/frames-adc.tar.gz) )
    debug "hf two: ${#hf[@]}: ${hf[*]}"
    [[ "${#hf[@]}" -eq 2 ]]
    debug "current with last three:"
    hf=$(historical_files --current -l 3 test-addnoise/frames-adc.tar.gz)
    debug "hf: $hf"
    hf=( $hf )
    debug "hf:: ${#hf[@]}: ${hf[*]}"
    [[ "${#hf[@]}" -eq 3 ]]
}
