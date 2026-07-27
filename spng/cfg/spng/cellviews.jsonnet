// This defines subgraphs configuring CellViews.
//
// Note, unlike the similar older CrossViews, this outputs one tensor per target
// view which has both MP2 and MP3 images in dimension dim=-3.


local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";

function(control={}, pg=real_pg)
{
    /// A subgraph with one input and one output port passing ITorchTensorSet
    /// though CellViews.  
    cellviews_tensorset(tpc,
                        out_views=[0,1],   // order and views to output
                        chunk_size=0,      // large detectors may want to set this
                        uvw_index=[0,1,2], // order of per-view Boolean tensors in set
                        extra_name="")::
        pg.pnode({
            type: 'SPNGCellViews',
            name: tpc.name + extra_name,
            data: {
                anode: wc.tn(tpc.anode),
                face_idents: tpc.faces,
                uvw_index: uvw_index,
                out_views: out_views,
                chunk_size: chunk_size,
            },
        }, nin=1, nout=1, uses=[tpc.anode]),
    
        

    /// A subgraph with three input ports and N output ports.  Input ports
    /// accept, in order, U, V and W ITorchTensor.  Output ports provide
    /// ITorchTensor for each target view in given order.  Large detectors may
    /// want to limit memory usage by specifying a chunk_size.
    cellviews_tensors(tpc, out_views=[0,1], chunk_size=0, extra_name="")::
        local this_name = tpc.name + '_cellviews' + extra_name;
        local packer = pg.pnode({
            type: 'SPNGTensorPacker',
            name: this_name,
            data: {
                multiplicity: 3
            } + control
        }, nin=3, nout=1);
        local cellviews = $.cellviews_tensorset(tpc, out_views, chunk_size, extra_name=extra_name);
        local nout = std.length(out_views);
        local unpacker = pg.pnode({
            type: 'SPNGTensorUnpacker',
            name: this_name,
            data: {
                selections: [{index: ind} for ind in wc.iota(nout)],
            } + control
        }, nin=1, nout=nout);
        pg.pipeline([packer, cellviews, unpacker]),

}

/*
        std::vector<int> face_idents = {0,1};

        /// The indices for the U, V and W tensors in the input tensor set.
        std::vector<int> uvw_index = {0,1,2};

        /// Enumerate the views (plane indices) for which output tensors will be
        /// produced.
        std::vector<int> out_views = {0,1};

        /// Required, name of IAnodePlane component corresponding to the channels.
        std::string anode ="";

        /// Name and order of the kinds of "cell views" to produce.  Each cell
        /// view maps to one index on the dim=-3 dimension of the output
        /// tensors.  
        std::vector<std::string> cell_views = {"mp2", "mp3"};

        /// If set, perform chunked processing
        int chunk_size = 0;

*/
