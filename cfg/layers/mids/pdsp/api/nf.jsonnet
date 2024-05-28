// The .nf sub-API method returns a subgraph for WCT noise filtering.

local low = import "../../../low.jsonnet";
local tn = low.wc.tn;
local get = low.wc.get;

local chndbs = {
    "actual": import "chndb-actual.jsonnet",
    "perfect": import "chndb-perfect.jsonnet",
};
local filters = import "noise-filters.jsonnet";
local frs = import "frs.jsonnet";

function(services, params, options) function(anode, name)
    local pars = std.mergePatch(params, get(options, "params", {}));

    local chndb = chndbs[pars.nf.chndb](anode, pars.nf.binning,
                                        frs(pars).nf, services.dft);

    // filter look up
    local flu = filters(services, pars, anode, chndb);
    local fls = pars.nf.filters;

    local fobj = {
        channel: [flu[nam] for nam in fls.channel],
        grouped: [flu[nam] for nam in fls.grouped],
        status:  [flu[nam] for nam in fls.status],
    };

    // the NF
    low.pg.pnode({
        type: 'OmnibusNoiseFilter',
        name: name,
        data: {

            // Nonzero forces the number of ticks in the waveform
            nticks: 0,

            // channel bin ranges are ignored
            // only when the channelmask is merged to `bad`
            maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
            channel_filters:        [tn(obj) for obj in fobj.channel],
            grouped_filters:        [tn(obj) for obj in fobj.grouped],
            channel_status_filters: [tn(obj) for obj in fobj.status],
            noisedb: tn(chndb),
            intraces: '',
            outtraces: '',
            // intraces: 'orig' + ident,
            // outtraces: 'raw' + ident,
        },
    }, nin=1, nout=1, uses = fobj.channel + fobj.grouped + fobj.status)
    
