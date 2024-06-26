local wc = import "wirecell.jsonnet";
local nominal = import "nominal.jsonnet";

{
    resample_500ns : nominal + {
        binning : {
            tick: 500 * wc.ns,
            nticks: 6000,
        },

        ductor: super.ductor {
            field_file: "uboone-100ns.json.bz2",
        },
    },
    resample_512ns : nominal + {
        binning : {
            tick: 512 * wc.ns,
            nticks: 6000 * 512 / 500,
        },

        ductor: super.ductor {
            field_file: "uboone-64ns.json.bz2",
        },

    },
}
