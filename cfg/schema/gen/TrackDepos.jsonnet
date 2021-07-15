/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This holds constructor functions for creating objects from the
 * schema: gen.TrackDepos
 * 
 * Prefer these functions over manual object construction to assure
 * greater validity.
 *
 */
local util_Base = import "schema/util/Base.jsonnet";

util_Base +
{
    gen: { TrackDepos: {
    // Construct Track (record)
    // 
    Track(obj=null, time=0.0, charge=-1.0, ray={}) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            time: $.util.Base.Time(time),
            charge: $.util.Base.Charge(charge),
            ray: $.util.Base.Ray(ray),
        },
    }.res,

    // Construct Tracks (sequence)
    // A sequence of tracks
    Tracks(val) :: assert(std.type(val)=='array'); [$.gen.TrackDepos.Track(v) for v in val],

    // Construct Config (record)
    // Configuration for TrackDepos component
    Config(obj=null, step_size=1.0, clight=1.0, group_time=-1.0, tracks=[]) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            step_size: $.util.Base.Distance(step_size),
            clight: $.util.Base.Normalized(clight),
            group_time: $.util.Base.Time(group_time),
            tracks: $.gen.TrackDepos.Tracks(tracks),
        },
    }.res,

    } }, // gen.TrackDepos
}
