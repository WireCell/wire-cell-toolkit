// The .nf sub-API method returns a function taking anode that when
// called gives subgraph for WCT noise filtering.

local low = import "../../../low.jsonnet";
local wc = low.wc;
local pg = low.pg;
local chndb = import "details/chndb.jsonnet";

function(services, params) function(anode)
    local name = "%d"%anode.data.ident;

    local chndbobj = chndb(services, params, anode);

    local bitshift = {
        type: "mbADCBitShift",
        name:name,
        data: {
            Number_of_ADC_bits: params.digi.resolution,
            Exam_number_of_ticks_test: 500,
            Threshold_sigma_test: 7.5,
            Threshold_fix: 0.8,
        },
    };

    local status = {
        type: "mbOneChannelStatus",
        name:name,
        data: {
            Threshold: 3.5,
            Window: 5,
            Nbins: 250,
            Cut: 14,
            anode: wc.tn(anode),
            dft: wc.tn(services.dft),
        },
        uses: [anode, services.dft],
    };

    local single = {
        type: "mbOneChannelNoise",
        name:name,
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(services.dft),
        },
        uses: [anode, services.dft, chndbobj],
    };
    
    local grouped = {
        type: "mbCoherentNoiseSub",
        name:name,
        data: {
            noisedb: wc.tn(chndbobj),
            anode: wc.tn(anode),
            dft: wc.tn(services.dft),
        },
        uses: [anode, services.dft, chndbobj],
    };

    local obnf = pg.pnode({
        type: "OmnibusNoiseFilter",
        name:name,
        data : {

            // This is the number of bins in various filters
            nsamples : params.nf.nsamples,

            maskmap: { chirp: "bad", noisy: "bad" },
            channel_filters: [
                wc.tn(bitshift),
                wc.tn(single)
            ],
            grouped_filters: [
                wc.tn(grouped),
            ],
            channel_status_filters: [
                wc.tn(status),
            ],
            noisedb: wc.tn(chndbobj),
            intraces: "orig",
            outtraces: "quiet",
        }
    }, uses=[chndbobj, anode, single, grouped, bitshift, status], nin=1, nout=1);

    local pmtfilter = pg.pnode({
        type: "OmnibusPMTNoiseFilter",
        name:name,
        data: {
            intraces: "quiet",
            outtraces: "raw",
            anode: wc.tn(anode),
        }
    }, nin=1, nout=1, uses=[anode]);
    
    pg.pipeline([obnf, pmtfilter], name=name)
