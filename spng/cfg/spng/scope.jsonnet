// A function to help describe a detector "scope"
//
// A scope is a detector name and an optional subset of TPCs.
//
// A scope is one part of a larger "system" object.  See system.jsonnet

local detconf = import "detconf.jsonnet";
local pg = import "pgraph.jsonnet";
local detector =  import "detector.jsonnet";

{
    bundle(detname, tpc_idents=[]):: {
        detname: detname,
        det: detconf[detname],
        tpcs: if std.length(tpc_idents) == 0
              then self.det.tpcs
              else [self.det.tpc[detector.tpc_name(ident)] for ident in tpc_idents],
    }
}

