/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This holds constructor functions for creating objects from the
 * schema: gen.AnodePlane
 * 
 * Prefer these functions over manual object construction to assure
 * greater validity.
 *
 */
local util_Base = import "schema/util/Base.jsonnet";

util_Base +
{
    gen: { AnodePlane: {
    // Construct Face (record)
    // Interesting drift positions on one anode face
    Face(obj=null, anode=0.0, response=0.0, cathode=0.0) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            anode: $.util.Base.Distance(anode),
            response: $.util.Base.Distance(response),
            cathode: $.util.Base.Distance(cathode),
        },
    }.res,

    // Construct Faces (sequence)
    // One or two faces
    Faces(val) :: assert(std.type(val)=='array'); [$.gen.AnodePlane.Face(v) for v in val],

    // Construct Config (record)
    // Configuration for AnodePlane
    Config(obj=null, ident=0, wire_schema="", nimpacts=10, faces=[]) :: {
        assert(std.setMember(std.type(obj), ["null", "object"])),
        res: if std.type(obj) == 'object' then obj else {
            ident: $.util.Base.Count(ident),
            wire_schema: $.util.Base.TypeName(wire_schema),
            nimpacts: $.util.Base.Count(nimpacts),
            faces: $.gen.AnodePlane.Faces(faces),
        },
    }.res,

    } }, // gen.AnodePlane
}
