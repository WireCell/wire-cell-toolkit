// see README.org

local wc = import "wirecell.jsonnet";

local base = import "chndb-base.jsonnet";
local perfect = import "chndb-perfect.jsonnet";
local rms_cuts = import "chndb-rms-cuts.jsonnet";

local frs = import "../frs.jsonnet";

function(services, params, anode)
    local fr = frs(params).nf;

    local types = {

        wct: function(epoch="before") {
            type: "OmniChannelNoiseDB",
            name: "ocndb%s"%epoch,
            data : {dft: wc.tn(services.dft)}+
                   if epoch == "perfect"
                   then perfect(params, anode, fr)
                   else base(params, anode, fr, rms_cuts[epoch]),
            uses: [anode, fr, services.dft],
        },

        wcls: function(epoch="before") {
            type: "wclsChannelNoiseDB",
            name: "wclscndb%s"%epoch,
            data : base(params, anode, fr, rms_cuts[epoch]) {
                misconfig_channel: {
                    policy: "replace",
                    from: {gain:  params.misconfigured.gain,
                        shaping: params.misconfigured.shaping},
                    to:   {gain: params.elec.gain,
                        shaping: params.elec.shaping},
                },
            },
            uses: [anode, fr],    // pnode extension
        },

        wcls_multi: function(name="") {
            local bef = $.wcls("before"),
            local aft = $.wcls("after"),
            type: "wclsMultiChannelNoiseDB",
            name: name,
            data: {
                rules: [
                    {
                        rule: "runbefore",
                        chndb: wc.tn(bef),
                        args: params.nf.run12boundary
                    },
                    {
                        rule: "runstarting",
                        chndb: wc.tn(aft),
                        args: params.nf.run12boundary,
                    },
                    // note, there might be a need to add a catchall if the
                    // above rules are changed to not cover all run numbers.
                ],
            },
            uses: [bef, aft],           // pnode extension
        },
    };
    if params.nf.chndb_epoch == "dynamic" || params.nf.chndb_type == "wcls_multi"
    then types.wcls_multi()
    else types[params.nf.chndb_type](params.nf.chndb_epoch)



