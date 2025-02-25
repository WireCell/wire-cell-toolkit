
local high = import "layers/high.jsonnet";
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(detector="pdsp")
    local params = high.params(detector);
    local mid = high.api(detector, params);
    # normally, pg.uses() but that takes a "graph" object not a bare uses list.
    local anodes = mid.anodes();
    local uses = wc.unique_list(pg.resolve_uses(anodes));
    local plugins = ["WireCellApps","WireCellGen"];
    local appcfg = {
        type: "AnodeDumper",
        data: {
            anodes: [wc.tn(anode) for anode in anodes],
        },
    };
    local cmdline = {
        type: "wire-cell",
        data: {
            plugins: plugins,
            apps: [appcfg.type]
        }
    };
    [cmdline] + uses + [appcfg]

