local wc = import "wirecell.jsonnet";
local p = import "pgrapher/experiment/pdsp/params.jsonnet";
local adc = p.adc;
local vrange = adc.fullscale[1] - adc.fullscale[0];
adc {
    voltmin: adc.fullscale[0]/wc.volt,
    voltmax: adc.fullscale[1]/wc.volt,
    voltrange: self.voltmax-self.voltmin,
    voltpeds: [b/wc.volt for b in adc.baselines],
    adccounts: std.pow(2,adc.resolution),
    adcpeds: [std.floor(self.adccounts*(vp-self.voltmin)/self.voltrange) for vp in self.voltpeds],
    voltperadc: self.voltrange/self.adccounts,
    voltageperadc: (adc.fullscale[1]-adc.fullscale[0])/self.adccounts,
    adcpervolt: self.adccounts/self.voltrange,
} 
  
