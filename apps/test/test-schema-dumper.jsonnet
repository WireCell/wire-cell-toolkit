// Simple test configuration for SchemaDumper
// This demonstrates basic usage of the schema-dumper function
//
// Run with:
//   WIRECELL_PATH=cfg:build:apps/test \
//   build/apps/wire-cell -p WireCellApps -p WireCellGen -p WireCellSigProc \
//   -c apps/test/test-schema-dumper.jsonnet
//
// Output includes factory schema with INode details and default configurations:
// {
//   "factories": {
//     "AddNoise": {
//       "classname": "AddNoise",
//       "concrete_type": "WireCell::Gen::IncoherentAddNoise",
//       "interfaces": ["WireCell::IConfigurable", "WireCell::IFrameFilter"],
//       "default_configuration": {  // from IConfigurable::default_configuration()
//         "dft": "FftwDFT",
//         "model": "",
//         "nsamples": 9600,
//         "replacement_percentage": 0.02,
//         "rng": "Random"
//       },
//       "node": {  // from INode methods (input_types(), output_types(), etc.)
//         "category": "functionNode",
//         "input_types": ["WireCell::IFrame"],
//         "output_types": ["WireCell::IFrame"],
//         "signature": "WireCell::IFrameFilter",
//         "concurrency": 1
//       }
//     },
//     "DumpFrames": {  // Example: INode but not IConfigurable
//       "classname": "DumpFrames",
//       "concrete_type": "WireCell::Gen::DumpFrames",
//       "interfaces": ["WireCell::IFrameSink"],
//       "node": {
//         "category": "sinkNode",
//         "input_types": ["WireCell::IFrame"],
//         "concurrency": 1
//       }
//     },
//     "AnodePlane": {  // Example: IConfigurable but not INode
//       "classname": "AnodePlane",
//       "concrete_type": "WireCell::Gen::AnodePlane",
//       "interfaces": ["WireCell::IConfigurable", "WireCell::IAnodePlane"],
//       "default_configuration": {
//         "faces": [null, null],
//         "ident": 0,
//         "nimpacts": 10,
//         "wire_schema": ""
//       }
//     }
//   },
//   "metadata": {
//     "generator": "WireCell::SchemaDumper",
//     "num_factories": 144,
//     "num_interfaces": 65
//   }
// }
//
// Statistics for typical plugin set:
//   - 144 total factories
//   - 133 with default_configuration (IConfigurable)
//   - 84 with node info (INode)
//   - 76 with both node and config
//   - 8 only node (INode but not IConfigurable)
//   - 57 only config (IConfigurable but not INode)
//
local dumper = import "schema-dumper.jsonnet";

// Dump all factories from common plugins to stdout
dumper()
