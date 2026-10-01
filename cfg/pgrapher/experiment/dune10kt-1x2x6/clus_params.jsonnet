// clus_params.jsonnet -- DUNE FD-HD 1x2x6 constants of the clustering and pattern-recognition
// jobs (clus.jsonnet, wct-clustering.jsonnet, pr.jsonnet, wct-pr-perevt.jsonnet).
//
// Promoted from wcp-porting-validation wcfm/wcfm_params.jsonnet (wcfm/docs/02, 25); the wcfm
// campaign file now takes these fields from here, so there is one copy.
//
// Geometry (verified against dune10kt-1x2x6-wires-larsoft-v1.json.bz2, wcfm/docs/02 sec 0):
// 12 APAs, all at x = 0, laid out 2 rows in y (even ident y < 0, odd ident y > 0) x 6
// columns in z (2306.4 mm long, 2323.9 mm pitch); BOTH faces live -- face 0 (+x wires,
// W channels 2080-2559) drifts to the cathode at +3629 mm, face 1 (-x, W 1600-2079) to
// -3629 mm; U (0-799) and V (800-1599) channels are wrapped over both faces.
// So x = 0 is an ANODE plane here, not a cathode (PDHD and SBND have the cathode at x = 0).
local wc = import 'wirecell.jsonnet';
local base_maker = import 'pgrapher/experiment/dune10kt-1x2x6/simparams.jsonnet';

{
    // The params object the jobs are built from (the same one the simulation uses).
    params: base_maker({}),

    // Readout / drift.  lar.drift_speed is the common base's 1.6 mm/us; the field file
    // (dune-garfield-1d565) says 1.565.  Kept at 1.6 everywhere so the simulation
    // (Drifter), BlobDepoFill and the clustering x agree with each other (wcfm/docs/02 sec 7).
    drift_speed: $.params.lar.drift_speed,
    tick: $.params.daq.tick,          // 0.5 us
    nticks: $.params.daq.nticks,      // 6000
    tick_span: 4,                     // MaskSlices span (imaging slices of 4 ticks)

    // Per-face drift-volume x extents for the clustering DetectorVolumes metadata (FV_x*):
    // from the W wire plane (|x| = 30.0155 mm) to the cathode face (apa_cpa - cpa_thick/2).
    fv_x: {
        wire: 30.0155 * wc.mm,
        cathode: 3.63075 * wc.m - 0.5 * 3.175 * wc.mm,   // 3629.1625 mm
    },
    // Overall active box (wires file extents; the 15 cm insets are the PDHD convention).
    fv_overall: {
        FV_xmin: -$.fv_x.cathode, FV_xmax: $.fv_x.cathode,
        FV_ymin: -6001.2 * wc.mm + 15 * wc.cm, FV_ymax: 6001.2 * wc.mm - 15 * wc.cm,
        FV_zmin: 0 * wc.mm + 15 * wc.cm, FV_zmax: 13925.9 * wc.mm - 15 * wc.cm,
    },
    // The active box itself (no insets): the PR job's fiducial box, from which the taggers'
    // margins are subtracted (pr.jsonnet).
    active_box: {
        xmin: -$.fv_x.cathode, xmax: $.fv_x.cathode,
        ymin: -6001.2 * wc.mm, ymax: 6001.2 * wc.mm,
        zmin: 0 * wc.mm, zmax: 13925.9 * wc.mm,
    },

    // Drift groups for clustering: all run anodes seen through face 0 (+x volume) and all
    // through face 1 (-x volume).  PDHD's ident%2 grouping is WRONG here (ident parity is
    // the y row, not the drift side).  `anodes` = the tools.anodes objects of the run.
    face_groups(anodes):: [
        { name: 'groupf0', face: 0, anodes: anodes },
        { name: 'groupf1', face: 1, anodes: anodes },
    ],

    // Bee display detector tag.  wire-cell-bee3 has a 'dune10kt-1x2x6' class (commit 1c89a0b)
    // but it is not deployed on the live server, where an unknown tag silently draws
    // MicroBooNE (wcfm/docs/17 sec 6); the PDHD tag only names the geometry Bee draws.
    bee_detector: 'protodunehd',
}
