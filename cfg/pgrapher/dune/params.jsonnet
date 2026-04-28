local common_params = import "pgrapher/common/params.jsonnet";
local wc = import "wirecell.jsonnet";

common_params {
    
    adc : super.adc{
        // resolution: std.extVar("Nbit"),
        resolution: 14,
    },
    
    elec: super.elec {
        // gain: std.extVar("elecGain")*wc.mV/wc.fC,
        gain: 7.8*wc.mV/wc.fC,
    }

}
