// Make an incoherent GroupNoiseModel "map file".
//
// Here we assume the test-noise-spectra.jsonnet with 10 groups
// (numbered 1 to 10) and PDSP APA0 channel numbering.  The particualr
// noise grouping is totally bogus for the real PDSP detector.
//
// We make each "plane" channel range use a different spectra.
//
// Note at the lowest noise, quantization error dominates and the
// modeled noise adding a relatively large component.  In this test,
// starting at about group 3, the modeled spectra are comparable to
// the input spectra in the peak but the modeled noise spectra will
// still have a fatter tail.
local step=256;
[
    // Note the "plane" attribute is not standard and shoud/must be
    // ignored.
    {
        plane: "p%d"%n,
        groupID: n,
        channels: std.range(step*(n-1), (step*n)-1),
    },
    for n in std.range(1,10)
]
