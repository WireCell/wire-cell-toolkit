
# The intention is to run this test in multiple releases and compare across releases.
@test "generate simple noise for comparison with older releases" {

    local ver="$(wire-cell --version)"
    local outfile="test-addnoise-${ver}.tar.bz2"

    run wire-cell -l test-addnoise.log -L debug \
        -A outfile="$outfile" -c gen/test/test-addnoise.jsonnet

    echo "$output"
    [[ "$status" -eq 0 ]]
    [[ -s "$outfile" ]]
}




