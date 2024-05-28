// The spdir and other tests make use of "ssss" (simple sim sig smeared splat something something)
// This helps generate the variant configs.
//
// Run wirecell-gen morse-* to find long/tran smearing numbers that match the
// extra spread the sigproc induces.  CAVEAT: these results are indirectly tied
// to the field response used in that process.  Defaults here are from PDSP.

function(nominal,
         long=[
             2.691862363980221,
             2.6750200122535057,
             2.7137567141154055
         ],
         tran=[
             0.7377218875719689,
             0.7157764520393882,
             0.13980698710556544
         ])
{
    local smeared = nominal + {
        splat: super.splat + {
            "smear_long": long,
            "smear_tran": tran,
        }
    },

    // Used for test-morse-pdsp
    morse_nominal: nominal,

    // Used for test-ssss-pdsp
    ssss_nominal: nominal,
    ssss_smeared: smeared,

    // Used for test/scripts/spdir 
    spdir_lofr: smeared {
        detector_data: super.detector_data {
            fields: "%s-fields-lo.json.bz2" % nominal.detector_name
        },
    },
    spdir_hifr: smeared {
        detector: super.detector_data {
            fields: "%s-fields-hi.json.bz2" % nominal.detector_name
        },
    },
}
