// This file provides ingredients for making IO objects.

local fileio = import "fileio.jsonnet";
local torchio = import "spng/torchio.jsonnet";
local pg = import "pgraph.jsonnet";

{
    /// WCT file I/O layers the organization of bits as:
    ///
    /// - schema describes how units of serialized data are identified and connected to WCT objects
    ///
    /// - format names of serialization mechanism
    ///
    /// - container names the file-at-rest technology
    ///
    /// There is some constraints and freedoms in forming supported
    /// (schema,format,container) triples.  For example, the WCT "frame" data
    /// model schema can be saved to a stream format which fills zip/npz or tar
    /// with or without compression or it can be saved to HDF which locks in the
    /// same format and container.
    ///

    /// Map common file extensions to the "format".  The formats are named:
    ///
    /// - sio :: these files contain a stream of named content using WCT's support for Boost Streams.
    ///
    /// - hio :: these files implement key/value random access HDF5 files.
    ///
    /// - pkl :: these files implement key/value random access Pickle files.
    ///
    /// Note, the file format does not necessarily imply a particular content
    /// schema.  Also, in the case of "sio", the "file format" refers to the
    /// stream format which is a container of files that can sink to many end
    /// container file formats.  For all formats
    format_by_ext(filename)::
        local choices = {
            // These are all implemented with WCT's Stream support.
            npz: "sio",
            tgz: "sio",
            tar: "sio",
            "tar.gz": "sio",
            "tar.bz2": "sio",
            zip: "sio",

            // These are implemented in specific sink / tap modules in the hio
            // subpackage.
            hdf: "hio",
            h5: "hio",
            hdf5: "hio",

            // These are specially implemented in SPNG and mostly for debugging.  Best not to propagate.
            pkl: "pkl",
            pickle: "pkl",
        };
        choices[$.file_ext(filename)],
        
    file_ext(filename)::
        local parts = std.split(filename, ".");
        parts[std.length(parts)-1],


    /// Return a sink of an IFrame directly to sio format.
    frame_to_sio(filename, tags=[], digitize=false, dense=null)::
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

    /// Return a sink of an IFrame directly to hio format.
    frame_to_hio(filename, digitize=false,  gzip=1, chunk=[256,256])::
        pg.pipeline([
            pg.pnode({
                type: 'HDF5FrameTap', // should have made it a sink....
                name: filename,
                data: {
                    filename: filename,
                    digitize: digitize,
                    gzip: gzip,
                    chunk: chunk
                }
            }, nin=1, nout=1),
            pg.pnode({
                type: 'DumpFrames',
                name: filename,
            }, nin=1, nout=0)]),

    /// Kitchen sink function that returns an IFrame sink based on file
    /// extension.  This requires providing the union of all possible config
    /// params not all of which may be used for the selected sink.
    frame_sink(filename,        // sio+hio
               digitize=false,  // sio+hio
               tags=[], dense=null, // sio
               gzip=1, chunk=[256,256] // hio
              )::
        local fmt = $.format_by_ext(filename);
        local choices = {
            sio: $.frame_to_sio(filename, tags=tags, digitize=digitize, dense=dense),
            hio: $.frame_to_hio(filename, digitize=digitize, gzip=gzip, chunk=chunk)
        };
        choices[fmt],

    // What is next is in terms of the tensor sets, ttensors are ITorchTensorSet
    // tensor sets and itensors are wct ITensorSet.  There are conversions
    // between them and there is conversion from IFrame to ITorchTensorSet.
    // Conversion to ITensorSet is needed for sio and hio but ITorchTensorSet
    // can go directly to pkl.

    // Each path to a final sink requires/accepts different config and we can
    // not decide that config policy here.  It's left to the caller to compose.

    /// Convert ITorchTensorSet to ITensorSet with optional inclusion and/or exclusion rules.
    ttensors_to_itensors(name, include_rules=[], exclude_rules=[], control={})::
        pg.pnode({
            type: 'SPNGTdmToTensorSet',
            name: name,
            data: {
                include_rules: include_rules,
                exclude_rules: exclude_rules,
            } + control,
        }, nin=1, nout=1),
    
    /// Return sink of ITensorSet to sio format.  See ttensors_to_itensors().
    itensors_to_sio(filename, prefix="")::
        pg.pipeline([
            pg.pnode({
                type: "TensorFileSink",
                name: filename,
                data: {
                    outname: filename,
                    prefix: prefix,
                },
            }, nin=1, nout=0)]),
                
    /// Return sink of ITensorSet to hio format.  See ttensors_to_itensors().
    itensors_to_hio(filename, datapath_pattern="tensorsets/{ident}", gzip=0, chunks=[])::
        pg.pnode({
            type: "HioTensorSink",
            name: filename,
            data: {
                filename: filename,
                gzip: gzip,
                chunks: chunks,
                datapath_pattern: datapath_pattern
            },
        }, nin=1, nout=0),

    /// Return sink of ITorchTensorSet directly to pkl.
    ttensors_to_pkl(filename)::
        pg.pnode({
            type: 'SPNGTensorSetPickleSink',
            name: filename, 
            data: {
                filename: filename,
            }
        }, nin=1, nout=0),


    /// Kitchen sink function that returns an ITorchTensorSet sink based on file
    /// extension.  This requires providing the union of all possible config
    /// params not all of which may be used for the selected sink.
    ttensors_sink(filename,     // hio+sio+pkl
                  // for hio + sio
                  include_rules=[], exclude_rules=[],
                  // for sio
                  prefix="",
                  // for hio
                  datapath_pattern="tensorsets/{ident}", gzip=0, chunks=[],
                  control={})::
        local fmt = $.format_by_ext(filename);
        local tti = $.ttensors_to_itensors(filename, include_rules=include_rules,
                                           exclude_rules=exclude_rules, control=control);
        local choices = {
            sio: pg.pipeline([tti, $.itensors_to_sio(filename, prefix=prefix)]),
            hio: pg.pipeline([tti, $.itensors_to_hio(filename, datapath_pattern=datapath_pattern,
                                                     gzip=gzip, chunks=chunks)]),
            pkl: $.ttensors_to_pkl(filename),
        };
        choices[fmt],

    /// older functions, should remove

    // Depos have essentially one source and sink which can handle different technology.
    depo_source:: fileio.depo_file_source,
    depo_sink:: fileio.depo_file_sink,

    // Serialize WCT ITensorSet (NOT ITorchTensorSet) to WCT "tensor data mode"
    // file whihc is .npy and .json files in .tar or .zip stream.
    tensor_source:: fileio.tensor_file_source,
    tensor_sink:: fileio.tensor_file_sink,


    // A dubious abstraction to return a sink function for a given file extension. 
    frame_sink_by_extension(ext)::
        local fmt = $.format_by_ext(ext);
        local choices = {
            sio:: $.frame_array_sink,
            hio:: $.frame_hdf_sink,
        };
        choices[fmt],

    /// Return a sync based on the file name extension.
    frame_array_any_sink(filename, digitize=false)::
        local parts = std.split(filename, ".");
        local ext = parts[std.length(parts)-1];
        $.frame_sink_by_extension(ext)(filename, digitize=digitize),

    // A "frame tensor file" is a "tensor file" in WCT frame tensor data model
    // (similar but different than SPNG tensor data model).  This serializes
    // .npy and .json files in .tar or .zip (.npz) streams.  It adds
    // Frame<-->ITensor serialization to basic tensor_{source,sink}
    frame_tensor_source:: fileio.frame_tensor_file_source,
    frame_tensor_sink:: fileio.frame_tensor_file_sink,

    // An older, simpler format for frames following WCT array data model.  This
    // serializes .npy files in .tar or .zip (.npz) streams.
    frame_array_source:: fileio.frame_file_source,
    frame_array_sink:: fileio.frame_file_sink,

    // A non-file sink.  It merely spews a bit of summary to the debug log.
    frame_null_sink(name, digitize=null):: pg.pnode({
        type: "DumpFrames",
        name: name,
        data: {},
    }, nin=1, nout=0),

    // This maps the WCT frame array data model to HDF.  There is currently only
    // a sink.  DNNROI training consumes files of this type.  Note, it does
    // substantial processing of the data prior to saving so does not exactly
    // represent IFrame.  It has several options.  Using default gzip=1 is
    // generally a good trade-off.  For ADCs it gives almost 2x smaller with
    // almost 2x writing time. For signals, it reduces size by 20x with less
    // than 2x writing time.  Using digitize=true for ADC is probably a good
    // idea.  Chunk sizes only matter for gzip != 0.  See
    // hio/docs/frame-npz-to-hio.org for one study on the matrix of speed and
    // compression in gzip vs chunk size.
    frame_hdf_sink(filename, digitize=false,  gzip=1, chunk=[256,256])::
        pg.pipeline([
            pg.pnode({
                type: 'HDF5FrameTap', // should have made it a sink....
                name: filename,
                data: {
                    filename: filename,
                    digitize: digitize,
                    gzip: gzip,
                    chunk: chunk
                }
            }, nin=1, nout=1),
            pg.pnode({
                type: 'DumpFrames',
                name: filename,
            }, nin=1, nout=0)]),


    // There is not really a choice for depos, this goes in all
    depos:: {
        depo_source:: $.depo_source,
        depo_sink:: $.depo_sink,
    },
        

    // Define some ready-to-use full "io objects"

    tensors:: $.depos {
        frame_source:: $.frame_tensor_source,
        frame_sink:: $.frame_tensor_sink,
    },

    arrays::  $.depos {
        frame_source:: $.frame_array_source,
        frame_sink:: $.frame_array_sink,
    },

}
