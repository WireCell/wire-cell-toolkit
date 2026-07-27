local wc = import "wirecell.jsonnet";
local real_pg = import 'pgraph.jsonnet';
local tpc_mod = import "spng/tpc.jsonnet";
local sg_js = import "spng/subgraphs.jsonnet";
local frame_js = import "spng/frame.jsonnet";
function(tpc, control, outname='test.pt', pg=real_pg) {

    local sg = sg_js(tpc, control.config, pg),
    local frame = frame_js(control),

    local tpc_nodes = tpc_mod(tpc, control, pg=pg),
    // local bypass = tpc_nodes.frame_bypass,
    local decon = tpc_nodes.response_decon,
    local wiener_filter = tpc_nodes.time_filter("wiener"),
    local unpack = tpc_nodes.frame_set_unpack,
    local repack = tpc_nodes.frame_set_repack,
    
    local repack1 = pg.pnode({
        type: 'SPNGTensorPacker',
        name: 'repack1',
        data: {
            multiplicity: 1,
        } + control
    }, nin=1, nout=1),

    // local from_tdm = tpc_nodes.frame_from_tdm,
    // local to_tdm = tpc_nodes.frame_to_tdm,
    // local to_tdm = sg.frame_to_tdm(extra_name="_TOTDM"),
    // local from_tdm = sg.tdm_to_frame(extra_name="_FROMTDM"),
    local to_tdm = frame.to_tdm(tpc),
    local from_tdm = frame.from_tdm(tpc, traces_tag="wiener"),
    local concat = pg.pnode({
            type: 'SPNGReduce',
            name: 'join',
            data: {
                multiplicity: 3,
                operation: "cat",
                tag: "asdf",
                dim: -2         // concatenate along channel dimension
            } + control
        }, nin=3, nout=1),
    local save = pg.pnode({
        type: 'SPNGTensorSetPickleSink',
        name: 'spng_save_tensor',
        data: {
            filename: outname
        }
    }, nin=1, nout=0),
    local fanout = pg.pnode({
        type: 'SPNGFanoutTensorSets',
        name: 'spng_fanout',
        data: {
            multiplicity: 2,
        }
    }, nin=1, nout=2),
    local fanin = pg.pnode({
        type: 'SPNGFaninTensorSets',
        name: 'spng_fanin',
        data: {
            multiplicity: 2,
        }
    }, nin=2, nout=1),

    // local rbl = [pg.pnode({
    //     type: 'SPNGRebaseliner',
    //     name: 'rebaseline%d'%i,
    //     data: {
    //         tag: "",
    //         datapath_format: "/traces/Rebaseliner/rebaseline%d"%i
    //     } + control
    // }, nin=1, nout=1) for i in std.range(0,2)],
    // local rbl_node = pg.intern(innodes=rbl, outnodes=rbl),

    // local decon_and_filter = pg.shuntlines([unpack, decon, wiener_filter, concat]),
    // local guts = pg.pipeline([decon_and_filter, repack1]),
    // local body = sg.wrap_bypass(guts),
    local decon_and_filter = pg.shuntlines([unpack, decon, wiener_filter, repack]),
    local body = sg.wrap_bypass(decon_and_filter),
    subgraph : pg.pipeline([to_tdm, body, from_tdm]),
    // subgraph : pg.intern(
    //     innodes=[to_tdm],
    //     outnodes=[from_tdm],
    //     centernodes=[unpack, decon, wiener_filter, repack, concat, repack1, save],
    //     edges=[
    //         pg.edge(to_tdm, fanout),
    //         pg.edge(fanout, unpack, 0),
    //         pg.edge(fanout, fanin, 1, 0),
    //         pg.edge(unpack, decon, 0, 0),
    //         pg.edge(unpack, decon, 1, 1),
    //         pg.edge(unpack, decon, 2, 2),
    //         pg.edge(unpack, decon, 3, 3),
    //         pg.edge(decon, wiener_filter, 0, 0),
    //         pg.edge(decon, wiener_filter, 1, 1),
    //         pg.edge(decon, wiener_filter, 2, 2),
    //         pg.edge(wiener_filter, concat, 0, 0),
    //         pg.edge(wiener_filter, concat, 1, 1),
    //         pg.edge(wiener_filter, concat, 2, 2),
    //         pg.edge(concat, repack1),
    //         pg.edge(repack1, fanin, 0, 1),
    //         pg.edge(fanin, from_tdm),


    //     ],
    // )

    // local test_shunt = pg.shuntlines([
    //     unpack, repack
    // ]),

    // subgraph : pg.intern(
    //     innodes=[to_tdm],
    //     outnodes=[from_tdm],
    //     centernodes=[bypass],
    // )
}