// The "system" object collects objects describing mutually orthogonal
// configuration which are also orthogonal to a particular job.
// - scope names detector and select TPCs
// - control give dft, rng, compute device, etc
// - io gives sources/sinks by data type making a choice of serialization
//
// Each system component has its own <name>.jsonnet file giving details.

local io = import "io.jsonnet";
local control = import "control.jsonnet";
local scope = import "scope.jsonnet";

{
    // Export Jsonnet modules to save user some lines.
    io: io, control: control, scope: scope,

    // Make a "system object".
    bundle(the_scope, the_control, the_io):: {
        scope: the_scope, control: the_control, io: the_io
    },

    // A vanilla / opinionated selection.  
    nominal(detname):: $.bundle(scope.bundle(detname), control.bundle(), io.arrays),
}


