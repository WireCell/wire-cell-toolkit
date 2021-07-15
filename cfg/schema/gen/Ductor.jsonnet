/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This holds constructor functions for creating objects from the
 * schema: gen.Ductor
 * 
 * Prefer these functions over manual object construction to assure
 * greater validity.
 *
 */
local util_Base = import "schema/util/Base.jsonnet";

util_Base +
{
    gen: { Ductor: {
    // Construct Config (record)
    // Ductor configuration
    Config(obj=null, nsigma=3, fluctuate=true, start_time=0.0, readout_time=5000000.0, tick=500.0, continuous=true, fixed=false, drift_speed=0.0, first_frame_number=0, anode="AnodePlane", rng="Random", pirs=[]) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            nsigma: $.util.Base.Count(nsigma),
            fluctuate: $.util.Base.Flag(fluctuate),
            start_time: $.util.Base.Time(start_time),
            readout_time: $.util.Base.Time(readout_time),
            tick: $.util.Base.Tick(tick),
            continuous: $.util.Base.Flag(continuous),
            fixed: $.util.Base.Flag(fixed),
            drift_speed: $.util.Base.Speed(drift_speed),
            first_frame_number: $.util.Base.Count(first_frame_number),
            anode: $.util.Base.TypeName(anode),
            rng: $.util.Base.TypeName(rng),
            pirs: $.util.Base.TypeNames(pirs),
        },
    }.res,

    } }, // gen.Ductor
}
