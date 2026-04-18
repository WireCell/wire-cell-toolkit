// Jsonnet function to configure wire-cell to dump factory schema information
//
// The SchemaDumper walks all registered interface types and collects information about
// registered factories by instantiating each one and querying its properties:
//
//   - classname: The registration key
//   - concrete_type: The actual C++ class name
//   - interfaces: List of all interfaces implemented
//
//   - default_configuration: (if factory is an IConfigurable)
//       The JSON/Jsonnet configuration object returned by default_configuration().
//       This shows all configuration parameters with their default values.
//
//   - node: (if factory is an INode) Node-specific information:
//       - category: sourceNode, sinkNode, functionNode, faninNode, fanoutNode,
//                   joinNode, splitNode, queuedoutNode, hydraNode
//       - input_types: List of C++ types for input ports (from input_types())
//       - output_types: List of C++ types for output ports (from output_types())
//       - signature: Signature string for the node
//       - concurrency: Number of concurrent instances allowed
//
// Note: INode and IConfigurable are independent - a factory may implement:
//   - Both (e.g., most filter nodes)
//   - Only INode (e.g., simple sinks like DumpFrames)
//   - Only IConfigurable (e.g., configuration-only components like AnodePlane)
//   - Neither (rare)
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
