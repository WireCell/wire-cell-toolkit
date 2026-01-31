// Simple test configuration for SchemaDumper
// This demonstrates basic usage of the schema-dumper function
//
// Run with:
//   wire-cell -c apps/test/test-schema-dumper.jsonnet
//
local dumper = import "schema-dumper.jsonnet";

// Dump all factories from common plugins to stdout
dumper()
