local wc = import "wirecell.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local sim = import "spng/sim.jsonnet";

local det = detconf.pdhd;
//det.tpcs[0].anode
local tpc = det.tpcs[0];
local anode = tpc.anode;
local face = if std.type(anode.data.faces[0]) == "null"
             then anode.data.faces[1]
             else anode.data.faces[0];

"%f cm" % [sim(det.tpcs[0], {}).response_distance / wc.cm]




