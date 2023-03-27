// This provides the variant "overlay" for *uboone* parameters.

// Do NOT place any "if" in this file!  If a new variant is needed,
// make an empty, new file in this directory and either import/inherit
// from this structure overiding only the new values.  Do not
// copy-paste.  Refactor this file if needed to achieve breif
// variants.
local wc = import "wirecell.jsonnet";

local nominal = import "nominal.jsonnet";

//
// gain_fudge_factor: 0.826
// sys_resp: false
//
// // Diffusion constants
// DiffusionLongitudinal: 6.4  // unit: cm^2/s, configured in the jsonnet file
// DiffusionTransverse:   9.8  // unit: cm^2/s, configured in the jsonnet file
//
// // Electron lifetime
// ElectronLifetime: 1000 // unit: ms, configured in the jsonnet file
//

nominal {
    ductor : super.ductor {

        // BIG FAT WARNING.  In "pure" WCT never never use std.extVar.
        // Until we can get TLAs working in the WCLS art tool we have
        // no great recourse but to put this mess here.  Thus the
        // "overlay" variant WILL NOT WORK outside of something that
        // calls Jsonnet with external variables.
        overlay: {
            // "MC_YZCorr_mcc9_v1_efieldcorr.root",
            filenameMC: std.extVar("YZCorrfilenameMC"),
            // ["hCorr0","hCorr1","hCorr2"],
            histnames: std.extVar("YZCorrhistnames"),
            // [0.0045, 0.0045, 0.0045],
            scaleDATA_perplane: std.extVar("scaleDATA_perplane"),
            // [0.0040, 0.0043, 0.0040],
            scaleMC_perplane: std.extVar("scaleMC_perplane"),
        }
    }, 

    misconfigured : super.misconfigured {
        channels: std.range(2016,2095) + std.range(2192,2303) + std.range(2352,2399),
        gain: 4.7*wc.mV/wc.fC,
        shaping: 1.1*wc.us,
    },
}
