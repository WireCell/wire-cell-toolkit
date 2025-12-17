// This provides funcitons at the full detector context.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local drift_mod = import "spng/drift.jsonnet";
local detsim_mod = import "spng/detsim.jsonnet";
local deposplat_mod = import "spng/deposplat.jsonnet";
local frame_mod = import "spng/frame.jsonnet";
local spng_mod = import "spng/spng.jsonnet";
local fans_mod = import "spng/fans.jsonnet";

function(det, control)
{
    // Intern these
    det: det,
    control: control,

    local fans = fans_mod(control),
    local frame = frame_mod(control),

    /// [1]IDepoSet->IDepoSet[1] node. Drifts the depos.
    drifter:
        drift_mod(det, control).drifter,


    // [1]IDepoSet->IFrame[ntpcs] node.  Applies parameterized model combining
    // simulation and signal processing to produce "true signal" waveforms.
    // There is one oport per TPC.
    splatter:
        deposplat_mod(det, control),


    // [1]IDepoSet->IFrame[ntpcs] node.  Applies detector response to depos to
    // produce ADC waveforms.  There is one oport per TPC.
    inducer:
        detsim_mod(det, control),


    /// Note, to fanout a 3 view node to both OSP and SPNG, use fans.fanout_select().

    // [ntpcs]IFrame->IFrame[ntpcs] node.  Input ADC frames to Original Signal
    // Processing / OmnibusSigProc to produce signals.
    osp:
        // FIXME: osp_subgraph is a temporary construct and will change.
        pg.crossline([tpc.osp_subgraphs.osp for tpc in det.tpcs]),

    // [ntpcs]IFrame->IFrame[ntpcs] node.  Input ADC frames to Signal
    // Processing, Next Generation subgraph to produce signals.
    spng:
        pg.crossline([spng_mod(tpc, control) for tpc in det.tpcs]),
    
    // [1]IDepoSet -> IFrame[ntpcs].  Depos input to drift and detsim to produce ADC frames.
    depos_to_adc:
        pg.pipeline([$.drifter, $.inducer]),

    // [1]IDepoSet -> IFrame[ntpcs].  Depos input to drift and splat to produce "true" signal frames.
    depos_to_splat:
        pg.pipeline([$.drifter, $.splatter]),
    
    // [1]IDepoSet -> IFrame[ntpcs].  Depos input to drift and sim and SP to produce "osp" signal frames.
    depos_to_osp:
        pg.shuntlines([$.drifter, $.inducer, $.osp]),

    // [1]IDepoSet->IFrame[3, ntpcs].  Produces an object with these attributes:
    //
    // - depo_sink :: A sink with one IDepoSet iport.
    // - splat_source :: A source of ntpcs of splat signals.
    // - osp_source :: A source of ntpcs of osp signals.
    // - spng_source :: A source of ntpcs of spng signals.
    //
    kitchen_sink:
        local drift = $.drifter;
        local splat = $.splatter;
        local sim = $.inducer;
        local depo_fan = pg.pnode({
            type:'DepoSetFanout',
            name: det.name + "_depo",
            data: {
                multiplicity:2
            }
        },nin=1, nout=2);

        local fanout_inode(view_index, M) = {
            type: "FrameFanout",
            name: det.name + '_v' + std.toString(view_index),
            data: {
                multiplicity: M
            }
        };

        local ntpcs = std.length(det.tpcs);
        local nadc_consumers = 2;
        local adc_fans = fans.fanout_cross_gen(ntpcs, nadc_consumers, fanout_inode);
        // a 3-in sink
        local adc_terminal = adc_fans[0];
        // 3-out sources
        local osp_adc = adc_fans[1][0];
        local spng_adc = adc_fans[1][1];

        local osp = $.osp;
        local spng = $.spng;


        // This holds the depo consumers (splat+sim).
        local depo_sink = pg.intern(innodes=[drift],
                                    centernodes=[depo_fan, splat, sim],
                                    edges=[
                                        pg.edge(drift, depo_fan),
                                        pg.edge(depo_fan, splat, 0, 0),
                                        pg.edge(depo_fan, sim, 1, 0)
                                    ] + [
                                        pg.edge(sim, adc_terminal, tpcid, tpcid)
                                        for tpcid in wc.iota(ntpcs)
                                    ]);
        local splat_source = pg.intern(oports=splat.oports);

        local osp_source = pg.intern(centernodes=[osp_adc],
                                     outnodes=[osp],
                                     edges=[
                                         pg.edge(osp_adc, osp, tpcid, tpcid)
                                         for tpcid in wc.iota(ntpcs)
                                     ]);
        local spng_source = pg.intern(centernodes=[spng_adc],
                                     outnodes=[spng],
                                     edges=[
                                         pg.edge(spng_adc, spng, tpcid, tpcid)
                                         for tpcid in wc.iota(ntpcs)
                                     ]);
        {
            depo_sink : depo_sink,
            splat_source : splat_source,
            osp_source : osp_source,
            spng_source : spng_source,

            test: {
                osp:osp,
                adc_terminal: adc_terminal,
                osp_adc: osp_adc,
                spng_adc: spng_adc,
                adc_fans: adc_fans,
            },
        }
                                     
                                   
        

}
