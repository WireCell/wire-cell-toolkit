// Jsonnet function to configure wire-cell to dump factory schema information
//
// Usage examples:
//   local dumper = import "schema-dumper.jsonnet";
//   dumper()                                           // dump to stdout with default plugins
//   dumper("schema.json")                             // dump to file
//   dumper("schema.json", ["WireCellGen", "WireCellSigProc"])  // with specific plugins
//
local wc = import "wirecell.jsonnet";

function(filename="/dev/stdout", plugins=["WireCellApps", "WireCellGen", "WireCellSigProc", "WireCellSio"])
    local appcfg = {
        type: "SchemaDumper",
        data: {
            filename: filename,
        },
    };
    local cmdline = {
        type: "wire-cell",
        data: {
            plugins: plugins,
            apps: [appcfg.type]
        }
    };
    [cmdline, appcfg]
