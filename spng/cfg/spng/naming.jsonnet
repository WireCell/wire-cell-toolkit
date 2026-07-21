local wc = import "wirecell.jsonnet";

local wpid_name(wpid) wc.unpack_wpid(wpid).name;

local anode_name(anode) "a" + std.toString(anode.data.ident);

local detector_name(detector) detector.name;


