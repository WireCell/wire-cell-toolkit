/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This holds constructor functions for creating objects from the
 * schema: WireCell.Sio.Cfg.NumpyFrameSaver
 * 
 * Prefer these functions over manual object construction to assure
 * greater validity.
 *
 */
local util_base = import "schema/util/base.jsonnet";

util_base +
{
    sio: { numpyframesaver: {
    // Construct Config (record)
    // Configuration for NumpyFrameSaver component
    Config(obj=null, digitize=false, baseline=0.0, scale=1.0, offset=0.0, frame_tags=[], filename="wct-frame.npz") :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            digitize: $.util.base.Flag(digitize),
            baseline: $.util.base.Charge(baseline),
            scale: $.util.base.Scaling(scale),
            offset: $.util.base.Charge(offset),
            frame_tags: $.util.base.Tags(frame_tags),
            filename: $.util.base.FilenameNPZ(filename),
        },
    }.res,


    } }, // sio.numpyframesaver
}
