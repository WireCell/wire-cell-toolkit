local high = import "layers/high.jsonnet";
local low = import "layers/low.jsonnet";
local wc = high.wc;
local pg = high.pg;

// Produce config sequence to test ShadowGhosting.
//
// These TLAs are accepted
//
// - detector :: canonical detector name
// 
// - spframes :: pattern matching frame files for each anode.
// 
function(detector="pdsp",
         spframes="frames-%(tier)s-%(anode)s.npz",
         clusters="clusters-%(tier)s-%(anode)s.npz")

    local mid = high.mid(detector);

    local anodes = mid.anodes();
    local nanodes = std.length(anodes);
    local anode_iota = std.range(0, nanodes-1);

    // make a ShadowGhosting config
    local sging(anode) = pg.pnode({
        type: 'ShadowGhosting',
        name: low.util.idents(anode),
    }, nin=1, nout=1);

    local apipes = [
        local anode = anodes[aid];
        pg.pipeline([
            high.fio.frame_file_source(std.format(spframes, {
                tier:"sig", anode:aid})),
            mid.img(anode),
            sging(anode),
            high.fio.cluster_file_sink(std.format(clusters, {
                tier:"deg", anode: aid})),
            
        ]) for aid in anode_iota];

    local graph = pg.components(apipes);
    high.main(graph, "TbbFlow")
