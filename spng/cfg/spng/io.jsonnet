// This file provides ingredients for making IO objects.
// 
// An IO object provides methods named like:
//
// <type>_{source,sink}
//
// An IO object is "opinionated" as to the file format.
//

local fileio = import "fileio.jsonnet";
local torchio = import "torchio.jsonnet";
local pg = import "pgraph.jsonnet";

{
    // first, give individual source/sinks

    // Depos have essentially one source and sink which can handle different technology.
    depo_source: fileio.depo_file_source,
    depo_sink: fileio.depo_file_sink,

    // Serialize WCT ITensorSet (NOT ITorchTensorSet) to WCT "tensor data mode"
    // file whihc is .npy and .json files in .tar or .zip stream.
    tensor_source: fileio.tensor_file_source,
    tensor_sink: fileio.tensor_file_sink,

    // Frames have several source/sink nodes to choose from.  Each take a
    // filename.  Sinks take a "digitize" argument.


    // A "frame tensor file" is a "tensor file" in WCT frame tensor data model
    // (similar but different than SPNG tensor data model).  This serializes
    // .npy and .json files in .tar or .zip (.npz) streams.  It adds
    // Frame<-->ITensor serialization to basic tensor_{source,sink}
    frame_tensor_source: fileio.frame_tensor_file_source,
    frame_tensor_sink: fileio.frame_tensor_file_sink,

    // An older, simpler format for frames following WCT array data model.  This
    // serializes .npy files in .tar or .zip (.npz) streams.
    frame_array_source: fileio.frame_file_source,
    frame_array_sink: fileio.frame_file_sink,

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
    depos: {
        depo_source: $.depo_source,
        depo_sink: $.depo_sink,
    },
        

    // Define some ready-to-use full "io objects"

    tensors: $.depos {
        frame_source: $.frame_tensor_source,
        frame_sink: $.frame_tensor_sink,
    },

    arrays:  $.depos {
        frame_source: $.frame_array_source,
        frame_sink: $.frame_array_sink,
    },




    // todo:
    // - [ ] clusters
    // - [ ] ROOT formats
}
