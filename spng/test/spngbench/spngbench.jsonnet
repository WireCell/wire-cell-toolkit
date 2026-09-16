// Factored SPNG-vs-OSP compute benchmark configuration.
//
// This builds ONE stage of a factored pipeline so each can be timed in
// isolation and so OSP and SPNG consume the *same* ADC input:
//
//   stage="sim"  : depos -> drift -> detsim            -> ADC frame file(s)
//   stage="osp"  : ADC frame file(s) -> OSP            -> signal frame file(s)
//   stage="spng" : ADC frame file(s) -> SPNG           -> signal frame file(s)
//
// The sim job is run once; its ADC frame file feeds both the OSP and SPNG jobs.
// Later the ADC file may be replaced by one made from real detector data.
//
// The det.jsonnet module exposes the factored pieces:
//   depos_to_adc = drift -> detsim   ([1]IDepoSet -> IFrame[ntpcs])
//   osp / spng   = the SP chains     ([ntpcs]IFrame -> IFrame[ntpcs])
//
// GPU sharding (SPNG only; OSP has a single DNN and uses `device` as-is) is
// expressed as a device MATRIX applied by rewriting data.device on per-APA SPNG
// node configs (whose names start "tpc<N>"):
//
//   gpu_scheme="none"         : every node on `device` (single device).
//   gpu_scheme="transverse"   : device = gpu[ tpc_index % ngpu ].  All nodes of
//       one APA share a GPU; APA pipelines are independent, so with one event in
//       flight only one runs at a time and the GPU never OOMs.  Up to (#APA) GPUs.
//   gpu_scheme="longitudinal" : device = gpu[ stage_index % ngpu ].  Graph stages
//       (decon/filter/dnn/roi, by node-name keyword) are spread across GPUs;
//       tensors cross GPUs at stage boundaries (ContextBase::to()).
//
// The matrix is (ngpu, gpu_scheme); extend the `stages` list below for finer
// longitudinal control.
//
// TLAs:
// @param input      stage=sim: depo file.  stage=osp/spng: ADC frame file (base).
// @param output     Output file (base).  Per-TPC files get a -tpc<N> suffix when
//                   ntpcs>1; ntpcs==1 uses the base name as-is.
// @param model_file roiuniter TorchScript (.ts) model (SPNG forward).
// @param stage      "sim" | "osp" | "spng".
// @param detname    Supported detector, default "pdhd".
// @param device     Base torch device: "cpu", "gpu", "gpu0", "gpu1", ...
// @param gpu_scheme "none" | "transverse" | "longitudinal".
// @param ngpu       Number of GPUs the scheme may use.
// @param engine     "Pgrapher" (serial) or "TbbFlow".
// @param wc_cores   TbbDataFlowGraph.max_threads (TbbFlow only; 0 = unlimited).
// @param verbosity  SPNG verbosity.
//
// Torch intra-op threads are set at run time via OMP_NUM_THREADS, not here.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local control_mod = import "spng/control.jsonnet";
local detconf = import "spng/detconf.jsonnet";
local det_mod = import "spng/det.jsonnet";
local roiuniter = import "spng/spng-roiuniter.jsonnet";

// ---------------------------------------------------------------------------
// Device matrix (inline so this file is self-contained; extend `stages`).
// ---------------------------------------------------------------------------
local devices = {
    // Ordered longitudinal stages; first matching keyword wins, stage index is
    // the list position.
    stages:: [
        { name: "decon",  keys: ["_group_"] },
        { name: "filter", keys: ["_gauss_", "_wiener", "_dnnroi_PRE", "_rebin", "_unbin", "_scale_"] },
        { name: "dnn",    keys: ["_dnnroi_fwd_", "_dnnroi_pre_", "_dnnroi_post_", "_dnnroi_roi_"] },
        { name: "roi",    keys: ["_applyroi_", "_cat", "_pack", "_groups_", "_TOTDM", "_FROMTDM", "_bypass"] },
    ],

    gpu(id):: "gpu" + std.toString(id),

    local digits = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"],
    local leading_int(s) =
        local acc = std.foldl(function(a, c)
                        if a.done then a
                        else if std.member(digits, c) then { done: false, s: a.s + c }
                        else { done: true, s: a.s }, std.stringChars(s), { done: false, s: "" });
        if acc.s == "" then null else std.parseInt(acc.s),

    // APA/TPC index encoded as the leading integer of a "tpc<N>..." node name.
    tpc_of(name)::
        if std.startsWith(name, "tpc")
        then leading_int(std.substr(name, 3, std.length(name) - 3))
        else null,

    // Longitudinal stage index by first matching keyword, or null.
    stage_of(name)::
        local hits = [i for i in std.range(0, std.length($.stages) - 1)
                      if std.length([k for k in $.stages[i].keys
                                     if std.length(std.findSubstr(k, name)) > 0]) > 0];
        if std.length(hits) > 0 then hits[0] else null,

    // assign(name) -> device string, or null to leave the node's base device.
    // Only per-APA SPNG nodes (tpc-named) are re-assigned; OSP nodes are left.
    assigner(scheme, ngpu):: function(name)
        local t = $.tpc_of(name);
        if t == null then null
        else if scheme == "transverse" then $.gpu(t % ngpu)
        else if scheme == "longitudinal" then
            local s = $.stage_of(name);
            if s == null then $.gpu(0) else $.gpu(s % ngpu)
        else null,

    // Rewrite data.device on any node config with one, per assign().
    reassign(cfg, assign):: [
        local obj = cfg[i];
        if std.isObject(obj) && std.objectHas(obj, "data") && std.isObject(obj.data)
           && std.objectHas(obj.data, "device") && std.objectHas(obj, "name")
           && assign(obj.name) != null
        then obj { data+: { device: assign(obj.name) } }
        else obj
        for i in std.range(0, std.length(cfg) - 1)
    ],

    // Split the shared torch forward service (SPNGTensorForwardTS, which holds
    // the model) into one clone per device used by the forward nodes referencing
    // it, and rewire each forward node to the clone on its own device.  Without
    // this, sharding puts (say) gpu1 inputs into a gpu0-resident model -> a
    // "tensors on cuda:0 and cuda:1" device mismatch.  A no-op when the forward
    // nodes all share one device (single GPU / no sharding).
    local is_fsvc(o) = std.isObject(o) && std.objectHas(o, "type")
                       && o.type == "SPNGTensorForwardTS",
    local is_fnode(o) = std.isObject(o) && std.objectHas(o, "type")
                        && o.type == "SPNGTensorForward",
    split_forward_services(cfg)::
        local fnodes = [o for o in cfg if is_fnode(o) && std.objectHas(o.data, "device")];
        local svcs = [o for o in cfg if is_fsvc(o)];
        local devs = std.set([o.data.device for o in fnodes]);
        if std.length(svcs) == 0 || std.length(devs) <= 1 then cfg
        else
            local sfx(d) = "@" + d;
            local clones = std.flattenArrays([
                [svc { name: svc.name + sfx(d), data+: { device: d } } for d in devs]
                for svc in svcs]);
            local rewritten = [
                local o = cfg[i];
                if is_fsvc(o) then null
                else if is_fnode(o) && std.objectHas(o.data, "device")
                        && std.objectHas(o.data, "forward")
                then o { data+: { forward: o.data.forward + sfx(o.data.device) } }
                else o
                for i in std.range(0, std.length(cfg) - 1)
            ];
            [x for x in rewritten if x != null] + clones,
};

function(input,
         model_file="",
         output="spngbench-out.npz",
         stage="spng",
         detname='pdhd',
         device='cpu',
         gpu_scheme='none',
         ngpu=1,
         napa=1,
         apa=-1,
         engine='Pgrapher',
         wc_cores=1,
         timeline='',
         verbosity=0)

    assert stage == "sim" || stage == "osp" || stage == "spng"
           : "spngbench: 'stage' must be sim|osp|spng, got " + stage;

    // Number of APAs (per-APA pipelines) the job covers, clamped to the
    // detector's physical APA count.  tpcids [0 .. napa-1] select that many
    // per-APA pipelines; more APAs means more independent pipelines in the graph
    // (and more wire-cell cores can be usefully applied).
    local phys_napa = std.length(detconf[detname].tpcs);
    local nap = std.max(1, std.min(wc.intify(napa), phys_napa));
    // apa>=0 selects a SINGLE physical APA (process-parallel axis: each process
    // runs one APA's pipeline).  apa<0 (default) runs the first `napa` APAs in
    // one job (the intra-process wc-core / GPU-shard axes).
    local apai = wc.intify(apa);
    local tpcids = if apai >= 0 then [std.min(apai, phys_napa - 1)]
                   else std.range(0, nap - 1);

    local osp_dump_prefix = std.strReplace(output, ".npz", "") + "_osp_dump";

    local controls = control_mod(device=device, verbosity=wc.intify(verbosity));
    // Pass `device` so OSP's DNN-ROI TorchService uses the same device as SPNG.
    // Without it detconf.get defaults the OSP TorchService to "cpu", so OSP's DNN
    // silently ran on CPU even for device=gpu* -- making OSP-vs-SPNG GPU
    // comparisons unfair (and OSP's VRAM zero).  For a single-GPU-shard base, OSP
    // is not sharded and uses this base device.
    local det = detconf.get(detname, tpcids, device=device, sp_dump_prefix=osp_dump_prefix);
    local ntpcs = std.length(det.tpcs);

    local spng_maker = roiuniter(model_file, do_transpose=false);
    local guts = det_mod(det, controls.config, spng_maker=spng_maker);

    // Per-TPC file name: base for a single TPC, else base with -tpc<N>.
    local ptname(base, t) =
        if ntpcs == 1 then base
        else std.strReplace(base, ".npz", "-tpc" + std.toString(t) + ".npz");

    local graph =
        if stage == "sim" then
            // depos -> drift -> detsim -> per-TPC ADC frame sink(s).
            local source = io.depo_source(input);
            local chain = guts.depos_to_adc;
            local sinks = [io.frame_array_sink(ptname(output, t)) for t in wc.iota(ntpcs)];
            if ntpcs == 1
            then pg.pipeline([source, chain, sinks[0]])
            else pg.intern(innodes=[pg.pipeline([source, chain])], outnodes=sinks,
                           edges=[pg.edge(chain, sinks[t], t, 0) for t in wc.iota(ntpcs)])
        else
            // ADC frame file(s) -> OSP|SPNG -> per-TPC signal sink(s).
            local sources = [io.frame_array_source(ptname(input, t)) for t in wc.iota(ntpcs)];
            local chain = if stage == "osp" then guts.osp else guts.spng;
            local sinks = [io.frame_array_sink(ptname(output, t)) for t in wc.iota(ntpcs)];
            if ntpcs == 1
            then pg.pipeline([sources[0], chain, sinks[0]])
            else pg.intern(innodes=sources, centernodes=[chain], outnodes=sinks,
                           edges=[pg.edge(sources[t], chain, 0, t) for t in wc.iota(ntpcs)]
                               + [pg.edge(chain, sinks[t], t, 0) for t in wc.iota(ntpcs)]);

    local tbb_dfp = {
        type: "TbbDataFlowGraph",
        name: "",
        data: { max_threads: wc.intify(wc_cores), summary: 1 }
              + (if timeline == '' then {} else { timeline: timeline }),
    };
    local extra_uses = if engine == "TbbFlow" then [tbb_dfp] else [];

    local base_cfg = pg.main(graph, app=engine,
                             plugins=["WireCellSpng", "WireCellSigProc", "WireCellGen", "WireCellPytorch"],
                             uses=controls.uses + extra_uses);

    // Apply GPU sharding by rewriting per-APA SPNG node devices, then split the
    // shared forward service so each APA's model lives on its inputs' GPU.
    if gpu_scheme == "none" then base_cfg
    else devices.split_forward_services(
        devices.reassign(base_cfg, devices.assigner(gpu_scheme, wc.intify(ngpu))))
