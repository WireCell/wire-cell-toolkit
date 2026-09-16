// Produce a WCT "depo file" of ionization depositions from a single muon
// tracked by edep-sim, as an alternative "probe" for the spngdir workflow.
//
// The muon starts at `pos` going in direction `dir` (both in global WCT
// coordinates -- pass the line-track meta's p0 and dir_glb so the muon matches
// the ideal line's geometry), at MIP kinetic `energy`.  The chain is:
//
//   ParticleGun -> ParticleTracking(edep-sim) -> Ionization -> TrackSegmentSampler -> DepoFileSink
//
// TLAs (all strings, parsed here):
//   output       : output depo file name (.npz / .tar / .tar.bz2).
//   gdml         : edep-sim detector GDML.  MUST tag the active LAr volume as a
//                  sensitive detector (see edep-sim's inputs/example.gdml) and
//                  be positioned to contain the track, or no segments result.
//   pos          : "x,y,z" start position in WCT units (mm), e.g. meta p0.
//   dir          : "x,y,z" direction (need not be normalized), e.g. meta dir_glb.
//   energy       : muon kinetic energy in MeV (default 1000).
//   count        : number of events (default 1).
//   physics_list : optional Geant4 reference physics list.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(output, gdml, pos="0,0,0", dir="0,0,1", energy="1000", count="1", physics_list="")

    local posv = [std.parseJson(x) for x in std.split(pos, ",")];
    local dirv = [std.parseJson(x) for x in std.split(dir, ",")];

    local gun = pg.pnode({
        type: "ParticleGun",
        data: {
            particle: "muon",
            energy: std.parseJson(energy) * wc.MeV,
            energy_type: "kinetic",
            position: posv,                       // WCT units (mm)
            direction: { dir: dirv, max_angle: 0.0 },
            count: std.parseJson(count),
            event0: 1,
        },
    }, nin=0, nout=1);

    local track = pg.pnode({
        type: "ParticleTracking",
        data: {
            gdml: gdml,
            physics_list: physics_list,
            segments: true,                       // only the ionization segments
        },
    }, nin=1, nout=1);

    local ion = pg.pnode({ type: "Ionization" }, nin=1, nout=1);

    local samp = pg.pnode({
        type: "TrackSegmentSampler",
        data: { ionization: "quanta", step_size: 1.0 * wc.mm },
    }, nin=1, nout=1);

    local sink = pg.pnode({
        type: "DepoFileSink",
        name: output,
        data: { outname: output },
    }, nin=1, nout=0);

    local graph = pg.pipeline([gun, track, ion, samp, sink]);

    pg.main(graph, app="Pgrapher",
            plugins=["WireCellGen", "WireCellEdep", "WireCellSio"])
