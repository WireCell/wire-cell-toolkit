// Simple test configuration for SchemaDumper
// This demonstrates basic usage of the schema-dumper function
//
// Run with:
//   WIRECELL_PATH=cfg:build:apps/test \
//   build/apps/wire-cell -p WireCellApps -p WireCellGen -p WireCellSigProc \
//   -c apps/test/test-schema-dumper.jsonnet
//
// Output includes factory schema with INode details:
// {
//   "factories": {
//     "FactoryName": {
//       "classname": "FactoryName",
//       "concrete_type": "Namespace::ClassName",
//       "interfaces": ["Interface1", "Interface2"],
//       "node": {  // only for INode factories
//         "category": "functionNode",
//         "input_types": ["WireCell::IFrame"],
//         "output_types": ["WireCell::IFrame"],
//         "signature": "WireCell::IFrameFilter",
//         "concurrency": 1
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
local dumper = import "schema-dumper.jsonnet";

// Dump all factories from common plugins to stdout
dumper()
