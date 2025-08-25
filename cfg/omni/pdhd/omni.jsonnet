// The omni object(s) for ProtoDUNE Horizontal-Drift (PDHD) detector (and variants)
//
// The contents in the directory holding this file, besides this file itself,
// are copied from:
//
// https://github.com/DUNE/dunereco/tree/develop/dunereco/DUNEWireCell/pdhd
//
// With this git hash and info from last commit:
// 741c6eb00f35bb6d6b60f3455ea2f139d094c5e5
// Mon Aug 18 13:12:30 2025 -0500
// v10_10_00d00

local wc = import 'wirecell.jsonnet';
local omniapi = import "../../omni.jsonnet";


// Currently, this grafts the "old style" configuration into the omni api.  As
// such, some API arguments may be ignored.

local tools_maker = import 'tools.jsonnet';
local params = import 'simparams.jsonnet';
local sim_maker = import 'sim.jsonnet';
local perfect = import 'chndb-base.jsonnet';
local nf_maker = import 'nf.jsonnet';
local sp_maker = import 'sp.jsonnet';

local chndb_base = import 'chndb-base.jsonnet';

local tools = tools_maker(params);
local sp = sp_maker(params, tools, { sparse: true });


// Basic PDHD omni.      
local omni = omniapi + {
    name: "pdhd",
    detname: "pdhd",
    lar: super.lar + params.lar,
    adc: super.adc + params.adc,
    volumes: params.det.volumes,
    binning: super.binning {
        default_time_binning: { spacing: 500*wc.ns, number: params.daq.nticks, start: 0.0*wc.ns },
        // fixme: all the stupid time offsets for sim/sp/splat
    },

    // FR is dependent on anode.
    responses(anode, kind, binning=$.binning.sim) :: {

        // Note, for now we do not distinguish based on "kind".
        local fr_file = if anode.data.ident == 0
                        then $.detdata.fields[1] // bad APA1/ident0
                        else $.detdata.fields[0], // nominal

        fr: [{
            type: 'FieldResponse',
            name: $.detname,
            data: { filename: fr_file } }],

        er: [{
            type: 'ColdElecResponse',
            name: $.detname,
            data: {
                start: binning.start,
                tick: binning.spacing,
                nticks: binning.number,
                shaping : params.elec.shaping,
                gain: params.elec.gain,
                postgain: params.elec.postgain,
            }}],

        local rcr = {
            type: 'RCResponse',
            name: $.detname,
            data: {
                width: params.rc_resp.width,
            }
        },
        // "RCRC" response
        rc: [rcr, rcr]
    },
    

    local make_chndb(anode, field, chndb_data = {}) = {
        type: 'OmniChannelNoiseDB',
        name: anode.name,
        data: chndb_base(params, anode, field, anode.data.ident) {
            dft:wc.tn($.dft)
        } + chndb_data,
        uses: [anode, field, $.dft],
    },

    //
    noisefilter(anode, binning=$.binning.sp, chndb_data={}, onf_data={}) ::
        local resp = $.responses(anode, "sp", binning);
        local chndb = make_chndb(anode, resp.fr[0], chndb_data);
        local nf = nf_maker(params, anode, chndb, anode.data.ident, anode.data.name);
        nf + { data: onf_data },


    // This disregards some arguments.
    sigproc(anode, binning=$.binning.sp, osp_data={}) ::
        sp.make_sigproc(anode),
        


};

// SPNG PDHD omni
local spng = import "spng.jsonnet";
local omni_spng = omni + {
    name: "pdhd-spng",
    detname: "pdhd",

    // This disregards some arguments.
    sigproc(anode, binning=$.binning.sp, osp_data={}) ::
        local resp = $.responses(anode, "sp", binning);
        spng.sigproc(anode, resp, $.adc)

        
        

};


[omni, omni_spng]
