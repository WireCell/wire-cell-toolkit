local g = import "pgraph.jsonnet";

{
    // a multi-layer fanout-pipelines-fanin structure
    // nnodes: number of nodes per layer
    // multi: multiplicity of nodes for a layer
    // fout is forward counting
    // fin is backward counting
    // e.g. for dune-vd 1x6x6, one can choose nnodes=[1,6], multi=[6,6]
    // for dune-vd 1x8x14, one can choose nnodes=[1,8,16], multi=[8,2,7]
    multifanpipe :: function( fout, pipelines, fin,
    fout_nnodes=[1,8,16], fout_multi=[8,2,7],
    fin_nnodes=[1,8,16], fin_multi=[8,2,7],
    name='multifanpipe', outtags=[], tag_rules=null ) {
        local fout_nlayers = std.length(fout_multi),
        assert fout_nlayers >= 2 : "fout_nlayers should be >= 2",
        local fin_nlayers = std.length(fin_multi),
        assert fin_nlayers >= 2 : "fin_nlayers should be >= 2",
        local npipe = std.length(pipelines),
        assert npipe == fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1] :
            "fout layout error npipe=%d, "%npipe + "fout=%d"%(fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1]),
        assert npipe == fin_nnodes[fin_nlayers-1]*fin_multi[fin_nlayers-1] :
            "fin layout error npipe=%d, "%npipe + "fin=%d"%(fin_nnodes[fin_nlayers-1]*fin_multi[fin_nlayers-1]),

        // function to create nodes for one layer
        local fout_layer(ilayer,nnodes,nmulti) = {
            ret : [
                g.pnode({
                    type: fout,
                    name: name+"_fout_%d"%ilayer + "_%d"%inode,
                    data: {
                        multiplicity: nmulti,
                        tag_rules: [],
                    }}, nin=1, nout=nmulti
                ) for inode in std.range(0,nnodes-1)],
        }.ret,
        // nodes for all layers
        local fout_layers = [
            fout_layer(ilayer,
            fout_nnodes[ilayer],
            fout_multi[ilayer])
            for ilayer in std.range(0,fout_nlayers-1)
        ],
        // make edges to make a combo node
        local fout_node = g.intern(
            innodes = fout_layers[0],
            centernodes = if fout_nlayers == 2 then [] else std.flattenArrays([fout_layers[i] for i in std.range(1,fout_nlayers-2)]),
            outnodes = fout_layers[fout_nlayers-1],
            edges = std.flattenArrays(
                [
                    [
                        g.edge(
                            fout_layers[ilayer-1][std.floor(inode/fout_multi[ilayer-1])],
                            fout_layers[ilayer][inode],
                            inode%fout_multi[ilayer-1],
                            0)
                    for inode in std.range(0,fout_nnodes[ilayer]-1)]
                for ilayer in std.range(1,fout_nlayers-1)])
        ),
        
        // similarly build the multi-layer fan in combo node
        // note the backward layer counting
        local fin_layer(ilayer,nnodes,nmulti) = {
            ret : [
                g.pnode({
                    type: fin,
                    name: name+"_fin_%d"%ilayer + "_%d"%inode,
                    data: {
                        multiplicity: nmulti,
                        tags: outtags,
                        tag_rules: [tag_rules for irule in std.range(0,nmulti-1)],
                    }}, nin=nmulti, nout=1
                ) for inode in std.range(0,nnodes-1)],
        }.ret,
        local fin_layers = [
            fin_layer(ilayer,
            fin_nnodes[ilayer],
            fin_multi[ilayer])
            for ilayer in std.range(0,fin_nlayers-1)
        ],
        local fin_node = g.intern(
            innodes = fin_layers[fin_nlayers-1],
            centernodes = if fin_nlayers == 2 then [] else std.flattenArrays([fin_layers[i] for i in std.range(1,fout_nlayers-2)]),
            outnodes = fin_layers[0],
            edges = std.flattenArrays(
                [
                    [
                        g.edge(
                            fin_layers[ilayer][inode],
                            fin_layers[ilayer-1][std.floor(inode/fin_multi[ilayer-1])],
                            0,
                            inode%fin_multi[ilayer-1])
                    for inode in std.range(0,fin_nnodes[ilayer]-1)]
                for ilayer in std.range(1,fin_nlayers-1)])
        ),

        // connect comb_fan_out-piples-combo_fan_in
        ret : g.intern(
            innodes = [fout_node],
            centernodes = pipelines,
            outnodes = [fin_node],
            edges = [g.edge(fout_node,pipelines[n],n,0) for n in std.range(0,npipe-1)] + 
            [g.edge(pipelines[n],fin_node,0,n) for n in std.range(0,npipe-1)],
        ),
    }.ret,

    // similar as multifanpipe but jusnt fan-out then pipelines with ending sinks
    multifanout :: function( fout, pipelines,
    fout_nnodes=[1,8,16], fout_multi=[8,2,7],
    name='multifanout', tag_rules=[] ) {
        local fout_nlayers = std.length(fout_multi),
        assert fout_nlayers >= 2 : "fout_nlayers should be >= 2",
        local npipe = std.length(pipelines),
        assert npipe == fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1] :
            "fout layout error npipe=%d, "%npipe + "fout=%d"%(fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1]),

        // function to create nodes for one layer
        local fout_layer(ilayer,nnodes,nmulti) = {
            ret : [
                g.pnode({
                    type: fout,
                    name: name+"_fout_%d"%ilayer + "_%d"%inode,
                    data: {
                        multiplicity: nmulti,
                        tag_rules: [],
                    }}, nin=1, nout=nmulti
                ) for inode in std.range(0,nnodes-1)],
        }.ret,
        // nodes for all layers
        local fout_layers = [
            fout_layer(ilayer,
            fout_nnodes[ilayer],
            fout_multi[ilayer])
            for ilayer in std.range(0,fout_nlayers-1)
        ],
        // make edges to make a combo node
        local fout_node = g.intern(
            innodes = fout_layers[0],
            centernodes = if fout_nlayers == 2 then [] else std.flattenArrays([fout_layers[i] for i in std.range(1,fout_nlayers-2)]),
            outnodes = fout_layers[fout_nlayers-1],
            edges = std.flattenArrays(
                [
                    [
                        g.edge(
                            fout_layers[ilayer-1][std.floor(inode/fout_multi[ilayer-1])],
                            fout_layers[ilayer][inode],
                            inode%fout_multi[ilayer-1],
                            0)
                    for inode in std.range(0,fout_nnodes[ilayer]-1)]
                for ilayer in std.range(1,fout_nlayers-1)])
        ),

        // connect comb_fan_out-piples
        ret : g.intern(
            innodes = [fout_node],
            centernodes = pipelines,
            outnodes = [],
            edges = [g.edge(fout_node,pipelines[n],n,0) for n in std.range(0,npipe-1)]
        ),
    }.ret,
}
