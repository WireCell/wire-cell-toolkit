#!/usr/bin/env bats

# Test for https://github.com/WireCell/wire-cell-toolkit/issues/200
#
# In part this issue invovles a problem building against WCT so we put
# this test outside of sub-packages.

@test "assure build config built" {
    eval $(../wcb dumpenv|grep =)
    [[ -f "$BuildConfig" ]]
}

@test "assure build config installed" {
    eval $(../wcb dumpenv|grep =)
    [[ -n "$PREFIX" ]]
    [[ -n "$BuildConfig" ]]

    echo "Expect build config file: $PREFIX/include/$BuildConfig"
    [[ -f "$PREFIX/include/$BuildConfig" ]]
}

@test "external build with config" {
    # assume we run from build/
    eval $(../wcb dumpenv|grep =)
    [[ -n "$PREFIX" ]]
    [[ -n "$BuildConfig" ]]
    [[ -n "$CXX" ]]

    here=$(pwd)
    tmp="$(mktemp -d ${BATS_TMPDIR}/test-issue200-build-with-config-XXXXX)"
    echo "working in: $tmp"
    cd $tmp/

    cat <<EOF > check.cxx
#include "$BuildConfig"
#include <iostream>
int main() {
    std::cout << "version " << WIRECELL_VERSION << std::endl;
    return 0;
}
EOF
    run $CXX -I $PREFIX/include -o check check.cxx
    echo $output
    [[ "$status" -eq 0 ]]

    run ./check
    echo "$output"
    [[ "$status" -eq 0 ]]

    [[ -n "$(echo $output|grep version)" ]]

    cd $here
    rm -rf $tmp
}
