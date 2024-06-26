#!/usr/bin/env bats

# This tests the "cfg/layers/" config against the original "cfg/pgrapher/"

bats_load_library wct-bats.sh

@test "test layer actual wcls sim nf sp" {
    cd_tmp file

    wcsonnet -o pgrapher-wct-sim-ideal-sn-nf-sp.json \
             pgrapher/experiment/uboone/wct-sim-ideal-sn-nf-sp.jsonnet 
    wcsonnet -o layers-wcls-sim-nf-sp.json layers/omnijob.jsonnet \
             -A tasks=sim,nf,sp \
             -A input=input.npz \
             -A output=output.npz \
             -A detector=uboone \
             -A variant=actual

}
