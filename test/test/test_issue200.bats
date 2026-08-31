#!/usr/bin/env bats

# Test for https://github.com/WireCell/wire-cell-toolkit/issues/200
#
# Here, we pretend we build against WCT and make sure the build config
# header is found.

bats_load_library wct-bats.sh

@test "assure build config built" {
    wcb_env
    echo  "$BuildConfig" 

    [[ -n "$BuildConfig" ]]
    path="$(topdir)/build/$BuildConfig"
    echo "Expect build config file built at: $path"
    [[ -f "$path" ]]
}

@test "assure build config installed" {
    wcb_env
    [[ -n "$PREFIX" ]]

    if [ ! -d $PREFIX/include/WireCellUtil ] ; then
        echo "Skip test for $BuildConfig because WCT is not installed"
        return
    fi

    echo "Expect build config file installed at: $PREFIX/include/$BuildConfig"
    [[ -f "$PREFIX/include/$BuildConfig" ]]

    [[ -n "$CXX" ]]

    cd_tmp

    local src=check-external-build-with-config

    cat <<EOF > ${src}.cxx
#include "$BuildConfig"
#include <iostream>
int main() {
    std::cout << "version " << WIRECELL_VERSION << std::endl;
    return 0;
}
EOF
    check $CXX -I $PREFIX/include -o ${src} ${src}.cxx

    check ./$src

    [[ -n "$(echo $output|grep version)" ]]


}
