/** Support for defining component configuration object schema.

Each configurable component shall define a schema for a record it
expects to receieve as configuration.

Each schema is defined in a namespace patterened after the package
name and the component class name.

*/

local moo = import "moo.jsonnet";

// the main output object:
function (pkg, comp) {

    // Set path convention incorporating package and component names.
    // This reflects in things outside of codegen control so do not
    // change it unless you really know what you are getting into.
    path: ["WireCell"+pkg, "Cfg", comp],

    // A schema factory used to make types "in" a package
    // namespace.  Use this in a "working object" context to
    // define needed final schema using local object keys to refer
    // to intermediate schema.
    schema :: moo.oschema.schema(self.path),
        
    // Call this on the "working object" to produce final schema
    // sequence.  Note, if any types outside the path are referred to
    // by types in hier, they must be prepended to the type sequence
    // returned by build() prior to final output.
    build(hier) :: moo.oschema.sort_select(hier,  self.path),
}


