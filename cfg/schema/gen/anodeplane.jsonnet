/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This holds constructor functions for creating objects from the
 * schema: WireCell.Gen.Cfg.AnodePlane
 * 
 * Prefer these functions over manual object construction to assure
 * greater validity.
 *
 */
local util_base = import "schema/util/base.jsonnet";

util_base +
{
    gen: { anodeplane: {
    // Construct Face (record)
    // Interesting drift positions on one anode face
    Face(obj=null, anode=0.0, response=0.0, cathode=0.0) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            anode: $.util.base.Distance(anode),
            response: $.util.base.Distance(response),
            cathode: $.util.base.Distance(cathode),
        },
    }.res,

    // Construct Faces (sequence)
    // One or two faces
    Faces(val) :: assert(std.type(val)=='array'); [$.gen.anodeplane.Face(v) for v in val],

    // Construct Config (record)
    // Configuration for AnodePlane
    Config(obj=null, ident=0, wire_schema="", nimpacts=10, faces=[]) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            ident: $.util.base.Count(ident),
            wire_schema: $.util.base.TypeName(wire_schema),
            nimpacts: $.util.base.Count(nimpacts),
            faces: $.gen.anodeplane.Faces(faces),
        },
    }.res,


    } }, // gen.anodeplane
}
