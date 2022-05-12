/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This holds constructor functions for creating objects from the
 * schema: WireCell.Gen.Cfg.TrackDepos
 * 
 * Prefer these functions over manual object construction to assure
 * greater validity.
 *
 */
local util_base = import "schema/util/base.jsonnet";

util_base +
{
    gen: { trackdepos: {
    // Construct Track (record)
    // 
    Track(obj=null, time=0.0, charge=-1.0, ray={}) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            time: $.util.base.Time(time),
            charge: $.util.base.Charge(charge),
            ray: $.util.base.Ray(ray),
        },
    }.res,

    // Construct Tracks (sequence)
    // A sequence of tracks
    Tracks(val) :: assert(std.type(val)=='array'); [$.gen.trackdepos.Track(v) for v in val],

    // Construct Config (record)
    // Configuration for TrackDepos component
    Config(obj=null, step_size=1.0, clight=1.0, group_time=-1.0, tracks=[]) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            step_size: $.util.base.Distance(step_size),
            clight: $.util.base.Normalized(clight),
            group_time: $.util.base.Time(group_time),
            tracks: $.gen.trackdepos.Tracks(tracks),
        },
    }.res,


    } }, // gen.trackdepos
}
