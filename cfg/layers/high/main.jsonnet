// This file is part of WCT "high layer" configuration structure.  It
// provides a "main()" function used to return the final WCT
// configuration sequence given a previously constructed graph.

local pg = import "pgraph.jsonnet";

local basic_plugins = [
    "WireCellSio", "WireCellAux",
    "WireCellGen", "WireCellSigProc", "WireCellImg", "WireCellClus",
    "WireCellApps"];

local app_plugins = {
    'TbbFlow': ["WireCellTbb"],
    'Pgrapher': ["WireCellPgraph"],
};

local findsubstrtype = function(objs, ss)
    [o for o in objs if std.length(std.findSubstr(ss, o.type)) > 0];

local hassubstr = function(objs, ss)
    std.length(findsubstrtype(objs, ss)) > 0;
          
// Guess plugin names based on object types.
local detect_plugins = function(uses) 
    std.set( std.prune([
        if hassubstr(uses, 'Torch') then 'WireCellPytorch',
        if hassubstr(uses, 'DNNROIFinding') then 'WireCellPytorch',
    ]));

function(graph, app='Pgrapher', extra_plugins = [])
    local uses = pg.uses(graph);
    local plugins = std.set(basic_plugins + extra_plugins + app_plugins[app] + detect_plugins(uses));
    local appcfg = {
        type: app,
        data: {
            edges: pg.edges(graph)
        },
    };
    local cmdline = {
        type: "wire-cell",
        data: {
            plugins: plugins,
            apps: [appcfg.type]
        }
    };
    [cmdline] + pg.uses(graph) + [appcfg]

