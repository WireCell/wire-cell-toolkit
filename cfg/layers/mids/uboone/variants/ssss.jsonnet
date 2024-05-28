local wc = import "wirecell.jsonnet";

local nominal = import "nominal.jsonnet";

local smeared = nominal + {
    splat: super.splat + {

        // Run wirecell-gen morse-* to find these numbers that match the extra
        // spread the sigproc induces.  CAVEAT: these results are inderectly
        // tied to the field response used.
        // CAVEAT2 and FIXME: these from PDSP
        "smear_long": [
            2.691862363980221,
            2.6750200122535057,
            2.7137567141154055
        ],
      "smear_tran": [
          0.7377218875719689,
          0.7157764520393882,
          0.13980698710556544
      ]        
    }
};


{
    // Used for test-morse-pdsp
    morse_nominal: nominal,

    // Used for test-ssss-pdsp
    ssss_nominal: nominal,
    ssss_smeared: smeared,

    // Used for test/scripts/spdir 
    spdir_lofr: smeared {
        detector: super.detector {
            fields: "uboone-fields-lo.json.bz2"
        },
    },
    spdir_hifr: smeared {
        detector: super.detector {
            fields: "uboone-fields-hi.json.bz2"
        },
    },
}
