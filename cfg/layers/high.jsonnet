// The "high level" layer returns an API object given a detector name
// and optional variant name.


// A lookup by detector and variant.  See gen-mids.sh.
local mids = import "mids.jsonnet";
local base = import "mids/base/api.jsonnet";
local svcs = import "high/svcs.jsonnet";

local midapi = import "midapi.jsonnet";

{
    // Forward the various "standard" wirecell utillities as a
    // convenience to end user.
    wc :: import "wirecell.jsonnet",
    pg :: import "pgraph.jsonnet",

    // Ways to move data between file and WCT.
    fio :: import "high/fileio.jsonnet",

    // This gives main(graph), which produces the final WCT
    // configuration sequence.
    main :: import "high/main.jsonnet",

    // The service pack factory function.  By default, the default
    // service pack is passed in to the mid factory function.  If a
    // non-default service pack is wanted, the caller should call
    // services() with some options and pass the result when creating
    // the mid.
    services :: svcs,

    // All known mid-layer detectors
    mids :: mids,

    // The parameter object factory.
    // 
    // The parameter object is selected from the detector's mid-level
    // variants.jsonnet structure based on the variant name.
    //
    // All three arguments are saved into the selected parameter object with a
    // key formed by appending "_name", eg "detector_name".
    //
    // The detector mid API may use "params.structure_name" to supply an API
    // that results in variant DFP graph structure independent from the chosen
    // parameter variant represented by "params.variant_name".
    //
    // The following variant and structure names SHALL be implemented according
    // to these guiding definitions:
    //
    // - nominal :: In name only.  This is for defining an idealized
    // configuration lacking variance.  It shall avoid introducing magic numbers
    // such as arbitrary time offsets, bad or per channel variance and shall
    // consider a uniform detector and resulting DFP graph structure.
    // 
    // - actual :: As close to the needs of realism as can be.  This
    // configuration shall endeavor to produce a configuration that model the
    // real world detector to the extent possible.  Strive to construct 
    // "actual" parameter variant and structure that derive from nominal.
    //
    // Other variant and structure names are allowed for special purposes.
    params :: function(detector, variant="nominal", structure=null)
        mids[detector].variants[variant] {
            detector_name: detector,
            variant_name: variant,
            structure_name: if std.type(structure) == "null" then variant else structure,
        },

    // Return a mid-level API implementation inside a midapi API.
    api :: function(detector, params, services=svcs(), options={})
        local sv = svcs();      //  fixme: need a way to choose eg GPU svcs
        //local base = midapi(sv, params, options=options);
        local def = base(sv, params, options=options);
        local det = mids[detector].api(sv, params, options=options);
        local imp = std.mergePatch(def, det);
        midapi(imp),
}


