// Produce config for WC 3D imaging subgraph for microboone

// This mostly transcribes:
// https://github.com/HaiwangYu/wcp-porting-img/blob/main/wct-celltree-img.jsonnet

// There are lots of "magic" numbers here that likely should be in a
// variant parameter pack.


local low = import "../../../low.jsonnet";
local pg = low.pg;
local wc = low.wc;
local high = import "layers/high.jsonnet";

// Return a subgraph that spans from frame to cluster.
// Shape of subgraph is controlled by params.img.slice_strategy
function(services, params) function(anode)
    local ident = low.util.idents(anode);
    local img = low.img(anode);


/*
    //-------------- haiwang's

local img = import "img.jsonnet";

local celltreesource = pg.pnode({
    type: "CelltreeSource",
    name: "celltreesource",
    data: {
        filename: "celltreeOVERLAY.root",
        EventNo: 6501,
        // in_branch_base_names: raw [default], calibGaussian, calibWiener
        in_branch_base_names: ["calibWiener", "calibGaussian"],
        out_trace_tags: ["wiener", "gauss"], // orig, gauss, wiener
        in_branch_thresholds: ["channelThreshold", ""],
    },
 }, nin=0, nout=1);

local dumpframes = pg.pnode({
        type: "DumpFrames",
        name: 'dumpframe',
    }, nin=1, nout=0);

local magdecon = pg.pnode({
      type: 'MagnifySink',
      name: 'magdecon',
      data: {
        output_filename: "mag.root",
        root_file_mode: 'UPDATE',
        frames: ['gauss', 'wiener', 'gauss_error'],
        cmmtree: [['bad','bad']],
        summaries: ['gauss', 'wiener', 'gauss_error'],
        trace_has_tag: true,
        anode: wc.tn(anodes[0]),
      },
    }, nin=1, nout=1, uses=[anodes[0]]);

local waveform_map = {
      type: 'WaveformMap',
      name: 'wfm',
      data: {
        filename: "microboone-charge-error.json.bz2",
      },
      uses: [],
    };

local charge_err = pg.pnode({
      type: 'ChargeErrorFrameEstimator',
      name: 'cefe',
      data: {
        intag: "gauss",
        outtag: 'gauss_error',
        anode: wc.tn(anodes[0]),
	rebin: 4,  // this number should be consistent with the waveform_map choice
	fudge_factors: [2.31, 2.31, 1.1],  // fudge factors for each plane [0,1,2]
	time_limits: [12, 800],  // the unit of this is in ticks
        errors: wc.tn(waveform_map),
      },
    }, nin=1, nout=1, uses=[waveform_map, anodes[0]]);

local cmm_mod = pg.pnode({
      type: 'CMMModifier',
      name: '',
      data: {
        cm_tag: "bad",
        trace_tag: "gauss",
        anode: wc.tn(anodes[0]),
	start: 0,   // start veto ...
	end: 9592, // end  of veto
	ncount_cont_ch: 2,
        cont_ch_llimit: [296, 2336+4800 ], // veto if continues bad channels
	cont_ch_hlimit: [671, 2463+4800 ],
	ncount_veto_ch: 1,
	veto_ch_llimit: [3684],  // direct veto these channels
	veto_ch_hlimit: [3699],
	dead_ch_ncount: 10,
	dead_ch_charge: 1000,
	ncount_dead_ch: 2,
	dead_ch_llimit: [2160, 2080], // veto according to the charge size for dead channels
	dead_ch_hlimit: [2176, 2096],
	ncount_org: 5,   // organize the dead channel ranges according to these boundaries 
	org_llimit: [0   , 1920, 3840, 5760, 7680], // must be ordered ...
	org_hlimit: [1919, 3839, 5759, 7679, 9592], // must be ordered ...
      },
    }, nin=1, nout=1, uses=[anodes[0]]);

local frame_quality_tagging = pg.pnode({
      type: 'FrameQualityTagging',
      name: '',
      data: {
        trace_tag: 'gauss',
        anode: wc.tn(anodes[0]),
	nrebin: 4, // rebin count ...
	length_cut: 3,
	time_cut: 3,
	ch_threshold: 100,
	n_cover_cut1: 12,
	n_fire_cut1: 14,
	n_cover_cut2: 6,
	n_fire_cut2: 6,
	fire_threshold: 0.22,
	n_cover_cut3: [1200, 1200, 1800 ],
	percent_threshold: [0.25, 0.25, 0.2 ],
	threshold1: [300, 300, 360 ],
	threshold2: [150, 150, 180 ],
	min_time: 3180,
	max_time: 7870,
	flag_corr: 1,
      },
    }, nin=1, nout=1, uses=[anodes[0]]);

local frame_masking = pg.pnode({
      type: 'FrameMasking',
      name: '',
      data: {
        cm_tag: "bad",
        trace_tags: ['gauss','wiener'],
        anode: wc.tn(anodes[0]),
      },
    }, nin=1, nout=1, uses=[anodes[0]]);

local anode = anodes[0];
// multi slicing includes 2-view tiling and dead tiling
local active_planes = [[0,1,2],[0,1],[1,2],[0,2],];
local masked_planes = [[],[2],[0],[1]];
// single, multi, active, masked
function(slicing_strategy = "single")

    local imgpipe = {
        single: pg.pipeline([
            img.slicing(anode, anode.name, "gauss", 109, active_planes=[0,1,2], masked_planes=[],dummy_planes=[]), // 109*22*4
            img.tiling(anode, anode.name),
            img.solving(anode, anode.name),
            // img.clustering(anode, anode.name),
            img.dump_new(anode, anode.name, params.lar.drift_speed),
        ], "img-" + anode.name),
        active: pg.pipeline([
            img.multi_active_slicing_tiling(anode, anode.name+"-ms-active", "gauss", 4),
            img.solving(anode, anode.name+"-ms-active"),
            img.dump(anode, anode.name+"-ms-active", params.lar.drift_speed)]),
        masked: pg.pipeline([
            // img.multi_masked_slicing_tiling(anode, anode.name+"-ms-masked", "gauss", 109),
            img.multi_masked_2view_slicing_tiling(anode, anode.name+"-ms-masked", "gauss", 109),
            img.clustering(anode, anode.name+"-ms-masked"),
            img.dump(anode, anode.name+"-ms-masked", params.lar.drift_speed)]),
        multi: {
            local active_fork = pg.pipeline([
                img.multi_active_slicing_tiling(anode, anode.name+"-ms-active", "gauss", 4),
                img.solving(anode, anode.name+"-ms-active"),
                img.dump(anode, anode.name+"-ms-active", params.lar.drift_speed),
            ]),
            local masked_fork = pg.pipeline([
                // img.multi_masked_slicing_tiling(anode, anode.name+"-ms-masked", "gauss", 109),
                img.multi_masked_2view_slicing_tiling(anode, anode.name+"-ms-masked", "gauss", 109),
                img.clustering(anode, anode.name+"-ms-masked"),
                img.dump(anode, anode.name+"-ms-masked", params.lar.drift_speed),
            ]),
            ret: pg.fan.fanout("FrameFanout",[active_fork,masked_fork], "fan_active_masked"),
        }.ret,
    };
    local graph = pg.pipeline([
        celltreesource,
        // frame_quality_tagging, // event level tagging
        // cmm_mod, // CMM modification
        // frame_masking, // apply CMM
        charge_err, // calculate charge error
        // magdecon, // magnify out
        // dumpframes,
        imgpipe[slicing_strategy],
    ], "main");

    high.main(graph, "TbbFlow")
    
*/

    
