// Configure wire-cell CLI to run a job that reads depos from file
// applies simulation, signal processing and imaging and writes files
// from the stages.  See depo-ssi-viz.sh for calling script that adds
// some more.

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

// Produces a top-level function.  Must give a detector type and may
// give a variant.  The outpat gives a pattern used to form output
// files.  A %(tier)s is interpolated by a tier name ("adc", "sig",
// "img") for each major link in the processing chain and a %(anode)s
// is interpolated by the anode plane ident.  Valid output extensions
// include .npz, .zip. .tar, .tar.gz, .tgz, etc.  If a indepos is
// given, read depos from there, else generate some internally.
function(detector, variant="nominal",
         indepos=null,
         outdepos="depos-drifted.npz",
         frames="frames-%(tier)s-%(anode)s.npz",
         clusters="clusters-%(tier)s-%(anode)s.npz")

    local mid = high.mid(detector, variant, options={sparse:false});

    local source = if std.type(indepos) == "null"
                  then mid.sim.track_depos()
                  else high.fio.depo_file_source(indepos);

    local drifter = mid.drifter();

    local anodes = mid.anodes();
    local nanodes = std.length(anodes);
    local anode_iota = std.range(0, nanodes-1);

    // Catch clusters for depo filling
    local catcher = function(aid)
        local anode = anodes[aid];
        local aname = anode.data.ident;
        local bsf = pg.pnode({
            type:'ClusterFanout',
            name: 'catch%d'%aname,
            data: {
                multiplicity: 2
            }}, nin=1, nout=2);
        local dbf = pg.pnode({
            type:'BlobDepoFill',
            name: aname,
            data: {
                // fixme, breaks detector independence
                speed: 1.56*wc.mm/wc.us, 
                // fixme, as does this
                time_offset: 314*wc.us,
            }}, nin=2, nout=1);
        local cfs = high.fio.cluster_file_sink(std.format(clusters, {
            tier:"tru", anode: anode.data.ident}));
        pg.intern([bsf, dbf], centernodes=[cfs],
                  edges=[pg.edge(bsf,dbf,1,0),
                         pg.edge(dbf,cfs,0,0)],
                  iports=[bsf.iports[0], dbf.iports[1]],
                  oports=[bsf.oports[0]],
                  name=aname);

    local catchers = [catcher(aid) for aid in anode_iota];

    local inject_catcher = function(img, cat, aid)
        local edge = {
            tail: {node: 'BlobGrouping:%d'%aid, port: 0},
            head: {node: 'ChargeSolving:%d'%aid, port:0}
        };
        pg.insert_node(img, edge, cat, cat, 0, 0);

    local apipes = [
        local anode = anodes[aid];
        pg.pipeline([
            mid.sim.signal(anode),
            mid.sim.noise(anode),
            mid.sim.digitizer(anode),

            high.fio.tap('FrameFanout',
                         high.fio.frame_file_sink(
                             std.format(frames, {
                                 tier:"adc", anode: anode.data.ident}),
                             tags=["orig%d" % anode.data.ident])),

            mid.sigproc.nf(anode),
            mid.sigproc.sp(anode),
            
            high.fio.tap('FrameFanout',
                         high.fio.frame_file_sink(
                             std.format(frames, {
                                 tier:"sig", anode: anode.data.ident}),
                             tags=["gauss%d" % anode.data.ident])),

            inject_catcher(mid.img(anode), catchers[aid], aid),

            high.fio.cluster_file_sink(std.format(clusters, {
                tier:"img", anode: anode.data.ident})),
        ]) for aid in anode_iota];

    local body = pg.fan.fanout('DepoSetFanout', apipes, "work");

    local mfan = pg.pnode({
        type:'DepoSetFanout',
        name:'play',
        data: {
            multiplicity: nanodes,
        }}, nin=1, nout=nanodes);
    local mind = pg.intern([mfan], centernodes=catchers,
                           edges=[high.pg.edge(mfan, catchers[aid], aid, 1),
                                  for aid in anode_iota]);
    local butt = high.fio.depo_file_sink(outdepos);
    local neck = pg.fan.fanout('DepoSetFanout', [body, mind, butt], name='drifted');

    local graph = pg.pipeline([source, drifter, neck], "main");
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

