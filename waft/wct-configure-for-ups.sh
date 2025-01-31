#!/bin/bash


usage () {
    cat <<EOF
Configure WCT source for building against a UPS products.

  wct-configure-for-ups.sh install_directory

This assumes the UPS environment for the dependent products is
already "setup".  

If install_directory is outside of UPS control then you can stop
reading now.

If the install_directory is the string "ups" then the source will be
configured to install into $WIRECELL_FQ_DIR.

That may point to a pre-existing "wirecell" UPS product.  If it does,
reusing it will install files on top of any pre-exising ones without
check so be careful.  

To make a fresh "wirecell" UPS produce area one can "declare" it like:

  $ ups declare wirecell <version> \
       -f \$(ups flavor) \
       -q e14:prof \
       -r wirecell/<version> \
       -z /path/to/install/products \
       -U ups  \
       -m wirecell.table

You'll have to provide the wirecell.table yourself, likely by copying
it from an existing "wirecell" UPS product.

Then, the calling environment can be munged like:

  $ setup wirecell <version> -q e14:prof

UPS is such a great and simple system! /s

EOF
    exit
}

install_dir="$1" ; shift
if [ "$install_dir" = "ups" ] ; then
    install_dir="$WIRECELL_FQ_DIR"
fi

# force to pick up GCC from PATH 
wct_cc=${CC:-gcc}
wct_cxx=${CXX:-g++}
wct_fort=${FORT:-gfortran}
export CXXFLAGS='-Wno-misleading-indentation -Wno-int-in-bool-context -Wvla'
env CC=$wct_cc CXX=$wct_cxx FC=wct_fort \
    ./wcb configure \
    --with-tbb=no \
    --with-jsoncpp="$JSONCPP_FQ_DIR" \
    --with-jsonnet="$GOJSONNET_FQ_DIR" \
    --with-spdlog="$SPDLOG_FQ_DIR" \
    --with-spdlog-lib="$SPDLOG_LIB" \
    --with-eigen-include="$EIGEN_DIR/include/eigen3" \
    --with-root="$ROOT_FQ_DIR" \
    --with-fftw="$FFTW_FQ_DIR" \
    --with-fftw-include="$FFTW_INC" \
    --with-fftw-lib="$FFTW_LIBRARY" \
    --boost-includes="$BOOST_INC" \
    --boost-libs="$BOOST_LIB" \
    --boost-mt \
    --prefix="$install_dir" \
    --build-debug="-O2 -ggdb3 -w" 

