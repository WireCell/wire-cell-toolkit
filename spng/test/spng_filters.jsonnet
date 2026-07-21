local wc = import 'wirecell.jsonnet';

function(device="cpu") {

    ROI_loose_lf : {
        type: "LfFilter",
        name: "ROI_loose_lf",
        data: {
            max_freq: 0.001,
            tau: 2.0e-06,
            use_negative_freqs: false,
        },
    },
    ROI_tight_lf : {
        type: "LfFilter",
        name: "ROI_tight_lf",
        data: {
            max_freq: 0.001,
            tau: 1.6e-05,
            use_negative_freqs: false,
        },
    },

    ROI_tighter_lf : {
        type: "LfFilter",
        name: "ROI_tighter_lf",
        data: {
            max_freq: 0.001,
            tau: 8.0e-050,
            use_negative_freqs: false,
        },
    },
    
    local gaus_filter = {
        type: 'HfFilter',
        name: 'Gaus_wide',
        data: {
            max_freq: 0.001,  // warning: units
            power: 2,
            flag: true,
            sigma: 0.00012,  // caller should provide
            use_negative_freqs: false,
            device: device,
        }
    },
    torch_gaus_filter : {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_gaus",
        uses: [gaus_filter],
        data: {
            spectra: [
                wc.tn(gaus_filter),
            ],
            device: device,
        },
    },

    local wiener_tight_u = {
        type: "HfFilter",
        name: "Wiener_tight_U",
        data: {
            flag: true,
            max_freq: 0.001,
            power: 6.55413,
            sigma: 0.000221933,
            use_negative_freqs: false,
        },

    },
    local wiener_tight_v = {
        type: "HfFilter",
        name: "Wiener_tight_V",
        data: {
            flag: true,
            max_freq: 0.001,
            power: 8.75998,
            sigma: 0.000222723,
            use_negative_freqs: false,
        },
    },

    local wiener_tight_w = {
        type: "HfFilter",
        name: "Wiener_tight_W",
        data: {
            flag: true,
            max_freq: 0.001,
            power: 3.47846,
            sigma: 0.000225567,
            use_negative_freqs: false,
        },
    },

    local wiener_wide_u = {
        type: "HfFilter",
        name: "Wiener_wide_U",
        data: {
            flag: true,
            max_freq: 0.001,
            power: 5.05429,
            sigma: 0.000186765,
            use_negative_freqs: false,
        },
    },

    local wiener_wide_v = {
        type: "HfFilter",
        name: "Wiener_wide_V",
        data: {
            flag: true,
            max_freq: 0.001,
            power: 5.77422,
            sigma: 0.0001936,
            use_negative_freqs: false,
        },
    },
    local wiener_wide_w = {
        type: "HfFilter",
        name: "Wiener_wide_W",
        data: {
            flag: true,
            max_freq: 0.001,
            power: 4.37928,
            sigma: 0.000175722,
            use_negative_freqs: false,
        },
    },
    local torch_wiener_wide_only_u = {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_wiener_wide_only_u",
        uses: [wiener_wide_u],
        data: {
            spectra: [
                wc.tn(wiener_wide_u),
            ],
            device: device,
        },
    },
    local torch_wiener_wide_only_v = {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_wiener_wide_only_v",
        uses: [wiener_wide_v],
        data: {
            spectra: [
                wc.tn(wiener_wide_v),
            ],
            device: device,
        },
    },
    local torch_wiener_wide_only_w = {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_wiener_wide_only_w",
        uses: [wiener_wide_w],
        data: {
            spectra: [
                wc.tn(wiener_wide_w),
            ],
            device: device,
        },
    },

    // wiener_wide_filters : [
    //     wiener_wide_u, wiener_wide_v, wiener_wide_w,
    // ],
    torch_wiener_wide_filters : [
        torch_wiener_wide_only_u, torch_wiener_wide_only_v, torch_wiener_wide_only_w,
    ],

    local torch_wiener_tight_only_u = {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_wiener_tight_only_u",
        uses: [wiener_tight_u],
        data: {
            spectra: [
                wc.tn(wiener_tight_u),
            ],
            device: device,
        },
    },
    local torch_wiener_tight_only_v = {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_wiener_tight_only_v",
        uses: [wiener_tight_v],
        data: {
            spectra: [
                wc.tn(wiener_tight_v),
            ],
            device: device,
        },
    },
    local torch_wiener_tight_only_w = {
        type: "Torch1DSpectrum",
        name: "torch_1dspec_wiener_tight_only_w",
        uses: [wiener_tight_w],
        data: {
            spectra: [
                wc.tn(wiener_tight_w),
            ],
            device: device,
        },
    },

    wiener_tight_filters : [
        wiener_tight_u, wiener_tight_v, wiener_tight_w,
    ],
    torch_wiener_tight_filters : [
        torch_wiener_tight_only_u, torch_wiener_tight_only_v, torch_wiener_tight_only_w,
    ],

    local wire_filters = [
        {
            type: 'HfFilter',
            name: 'Wire_ind',
            data: {
                max_freq: 1,  // warning: units
                power: 2,
                flag: false,
                sigma: 1.0 / wc.sqrtpi * 0.75,  // caller should provide
            }
        },
        {
            type: 'HfFilter',
            name: 'Wire_col',
            data: {
                max_freq: 1,  // warning: units
                power: 2,
                flag: false,
                sigma: 1.0 / wc.sqrtpi * 10.0,  // caller should provide
            }
        },
    ],

    torch_wire_filters : [
        {
            type: "Torch1DSpectrum",
            name: "torch_1dspec_ind",
            uses: [wire_filters[0]],
            data: {
                spectra: [
                    wc.tn(wire_filters[0]),
                ],
                device: device,
            },
        },
        {
            type: "Torch1DSpectrum",
            name: "torch_1dspec_col",
            uses: [wire_filters[1]],
            data: {
                spectra: [
                    wc.tn(wire_filters[1]),
                ],
                device: device,
            },
        },
    ],

}
