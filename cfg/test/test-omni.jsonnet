local omni = import "omni.jsonnet";
local omnimap = import "omni/all.jsonnet";



local getit(detname) =
    if detname == "base"
    then omni
    else omnimap[detname];

local testit(detname) = {
    local it = getit(detname),
    local anode=it.anodes[0],
    anode: anode,
    adc: it.adc,
    drift: it.drift(anode),
    splat: it.splat(anode),
    signal: it.signal(anode),
    noise: it.noise(anode),
    digitize: it.digitize(anode),
    nf: it.noisefilter(anode),
    sp: it.sigproc(anode),
};
    

{[name]:testit(name) for name in std.objectFields(omnimap)}

    
