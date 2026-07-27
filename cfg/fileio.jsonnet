/** General, detector file I/O methods to produce configuration objects.

A source or a sink takes a filename as first arg and this is used to provide a
name for the component instance.

The functions are named with a cross product of these labels: (type) x (format)
* ("sink" "source").  Where the type is "depo", "frame", "cluster" etc.  Format
is either a "direct" conversion format and has no label or is the label "tensor"
which indicates the use of the unified WCT tensor data model format.  The cross
product is "sparse" in that not all possible terms are implemented.

This file is taken from the experimental "layers" config.

See also omni.jsonnet.

*/

local pg = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";

{
    //
    // General tensor data model I/O nodes
    // 

    tensor_file_sink(filename, prefix="") ::
        pg.pnode({
            type: "TensorFileSink",
            name: filename,
            data: {
                outname: filename,
                prefix: prefix,
            },
        }, nin=1, nout=0),

    tensor_file_source(filename, prefix="") ::
        pg.pnode({
            type: "TensorFileSource",
            name: filename,
            data: {
                inname: filename,
                prefix: prefix,
            },
        }, nin=0, nout=1),

    ///
    /// Composite subgraph for tensor data model I/O with conversion
    /// between tensor and data product.
    ///

    ///
    ///  frame
    ///

    /// Sink a stream of frames to a file in WCT "frame tensor file
    /// format".  See aux/docs/frame-files.org for valid filename
    /// types and what the other options mean.  Prefer "frame tensor
    /// files" over the older, less capable "frame files".
    frame_tensor_file_sink(filename, prefix="",
                           mode="sparse", digitize=false) :: 
        pg.pipeline([
            pg.pnode({
                type: "FrameTensor",
                name: filename,
                data: {
                    digitize: digitize,
                    mode: mode,
                },
            }, nin=1, nout=1),
            pg.pnode({
                type: "TensorFileSink",
                name: filename,
                data: {
                    outname: filename,
                    prefix: prefix,
                },
            }, nin=1, nout=0)],
                   name=filename),
        
    /// Source a stream of frames from a file in WCT "frame tensor
    /// file format".
    frame_tensor_file_source(filename, prefix="") ::
        pg.pipeline([
            pg.pnode({
                type: "TensorFileSource",
                name: filename,
                data: {
                    inname: filename,
                    prefix: prefix,
                },
            }, nin=0, nout=1),
            pg.pnode({
                type: "TensorFrame",
                name: filename,
            }, nin=1, nout=1)],
                    name=filename),

    ///
    ///   cluster. See aux/docs/ClusterArrays.org
    ///

    /// Sink a stream of clusters to a file in WCT "cluster tensor
    /// file schema".  The datapath is used to "place" the main
    /// "cluster" tensor and will be interpolated with the cluster
    /// ident number.  The prefix is arbitrary.
    cluster_tensor_file_sink(filename, prefix="", datapath="clusters/%d") ::

        pg.pipeline([
            pg.pnode({
                type: "ClusterTensor",
                name: filename,
                data: {
                    datapath: datapath,
                },
            }, nin=1, nout=1),
            pg.pnode({
                type: "TensorFileSink",
                name: filename,
                data: {
                    outname: filename,
                    prefix: prefix,
                },
            }, nin=1, nout=0)], name=filename),
        
    /// Source a stream of clusters from a file in WCT "cluster tensor
    /// file schema".  An array of anode objects are required to
    /// provide anode faces matching idents referenced in the data.
    /// The datapath is a regular expression used to locate the
    /// "cluster" tensor.  The prefix must match any used to write the
    /// file.
    cluster_tensor_file_source(filename, anodes, prefix="", datapath="clusters/[[:digit:]]+$") ::
        pg.pipeline([
            pg.pnode({
                type: "TensorFileSource",
                name: filename,
                data: {
                    inname: filename,
                    prefix: prefix,
                },
            }, nin=0, nout=1),
            pg.pnode({
                type: "TensorCluster",
                name: filename,
                data: {
                    dpre: datapath,
                    anodes: [wc.tn(a) for a in anodes]
                }
            }, nin=1, nout=1, uses=anodes)], name=filename),


    // Direct converters


    // A source of depos from file. Several file formats supported.
    depo_file_source(filename, scale=1.0) ::
        pg.pnode({
            type: 'DepoFileSource',
            name: filename,
            data: { inname: filename, scale: scale }
        }, nin=0, nout=1),

    // A sink of depos.  Several file formats supported
    depo_file_sink(filename) ::
        pg.pnode({
            type: 'DepoFileSink',
            name: filename,
            data: { outname: filename }
        }, nin=1, nout=0),


    // Depreciated: prefer frame_tensor_file_sink().  Sink a stream
    // of frames to a file in WCT "frame file format".
    frame_file_sink(filename, tags=[],
                    digitize=false, dense=null) :: 
        pg.pnode({
            type: "FrameFileSink",
            name: filename,
            data: {
                outname: filename,
                tags: tags,
                digitize: digitize,
                dense: dense,
            },
        }, nin=1, nout=0),
        
    // Depreciated: prefer frame_tensor_file_source().  Source a
    // stream of frames from a file in WCT "frame file format".
    frame_file_source(filename, tags=[]) ::
        pg.pnode({
            type: "FrameFileSource",
            name: filename,
            data: {
                inname: filename,
                tags: tags,
            },
        }, nin=0, nout=1),


    // A cluster file is a tar stream, with optional compression,
    // holding json, dot or numpy format files.
    cluster_file_sink :: function(filename, format="json", prefix="cluster")
        pg.pnode({
            type: 'ClusterFileSink',
            name: filename,
            data: {
                outname: filename,
                format: format,
                prefix: prefix,
            },
        }, nin=1, nout=0),
    cluster_file_source :: function(filename, anodes, prefix="cluster")
        pg.pnode({
            type: 'ClusterFileSource',
            name: filename,
            data: {
                inname: filename,
                prefix: prefix,
                anodes: [wc.tn(a) for a in anodes],
            },
        }, nin=0, nout=1, uses=anodes),

    // A somewhat slippery format.
    celltree_file_source :: function(filename, recid, 
                                     branches = ["calibWiener", "calibGaussian"],
                                     frame_tags=["gauss"],
                                     trace_tags = ["wiener", "gauss"],
                                     extra_params = {})
        pg.pnode({
            type: "CelltreeSource",
            name: filename,
            data: {
                filename: filename,
                EventNo: std.toString(recid), // note, must NOT be a number but a string....
                frames: frame_tags,
                "in_branch_base_names": branches,
                "out_trace_tags": trace_tags,
            } + extra_params,
        }, nin=0, nout=1),


    // Adapt a sink to "tap" pattern via a 2-way fanout.  First output
    // port forwards input, second output port forwards to the sink.
    // The type should be a fanout type such as FrameFanout, etc.
    tap :: function(type, sink, tag_rules=[], name=null)
        local fanout = pg.pnode({
            type: type,
            name: pg.prune_array([name, sink.name, ""])[0],
            data: {
                multiplicity: 2,
                tag_rules: tag_rules,
            },
        }, nin=1, nout=2);
        pg.intern(innodes=[fanout],
                  outnodes=[fanout],
                  centernodes=[sink],
                  edges=[pg.edge(fanout, sink, 1, 0)],
                  name=name),

}
