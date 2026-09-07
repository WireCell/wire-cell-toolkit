// PDHD space-charge (curved) fiducial volume.  Counterpart of
// protodunevd/curved_fiducial.jsonnet (doc pdvd/41 sec 9, doc pdvd/43); the
// measurement behind PDHD's numbers is doc pdhd/09.
//
// The counterpart of the MicroBooNE prototype's ToyFiducial polygons
// (prototype_base/wire-cell/pid/src/ToyFiducial.cxx:60-135, 204-215): in each of the
// X-Y and X-Z planes the apparent detector boundary is not a box but a surface that is
// at the nominal wall near the anode and moves inward toward the cathode, and the two
// planes are tested INDEPENDENTLY and AND-ed (uBooNE: pnpoly(XY) && pnpoly(XZ), here
// CompositeFiducial{logic:'and'} of two PolyFiducial).
//
// TWO DELIBERATE DIFFERENCES FROM THE PDVD FORK
//
// 1. PDHD's y walls are NOT symmetric about zero.  PDVD's active volume is y = +-336.39,
//    so its side() calls pass -YW and +YW; PDHD's is y = 7.61 .. 606.0 (pr.jsonnet:1259-1260),
//    so the low and high walls are independent constants YLO and YHI.  Everything else in
//    the polygon assembly is structurally identical, because PDHD is also two drift volumes
//    mirrored about a cathode at x = 0 with anodes at +-XW.
//
// 2. There is no d50 default.  PDVD's function defaults are its doc-41 half-density
//    trapezoids; doc pdvd/41 sec 13.3 then WITHDREW that surface as a tagger boundary --
//    a fiducial for an ENDPOINT test must come from the endpoint distribution's tail, not
//    a charge-density median, and the d50 arm called 22.7 % of exit ends contained against
//    7.5 % for the flat box.  PDHD therefore never measured a d50 surface, and the defaults
//    below are dc = 0 on every wall: with no `profile` this file reproduces the NOMINAL
//    walls, i.e. pdhd_pr_fv's box.  A real surface only ever comes from `profile`.
//
// Each wall of each drift volume is one trapezoid in |x| (cm):
//
//     inset(|x|) = dc                          for |x| <= x1     (the cathode plateau)
//                = dc * (x2 - |x|)/(x2 - x1)   for x1 < |x| < x2 (the ramp)
//                = 0                           for |x| >= x2     (at the nominal wall)
//
// x2 <= the anode face is REQUIRED: the displacement a drifting charge accumulates is
// zero at the anode by construction.
//
// CUSHION.  These vertices are the calibrated surface itself, cushion 0.  MicroBooNE
// applies its cushion two different ways and this file supports both:
//   - the containment / physics fiducial insets the polygon VERTICES by
//     boundary_dis_cut = 3 cm (4 cm on z) -- pass cushion_y / cushion_z / cushion_x here;
//   - the cosmic tagger keeps the uncushioned surface and probes it with a tolerance
//     BAND (Cosmic_tagger.h:34-36) -- in WCT that is the taggers' fv_tolerance, which
//     offsets the POINT instead of the vertex (FiducialUtils.cxx:79-119).
// PDHD carries the SAME flat 15 cm space-charge allowance PDVD does, and this surface
// REPLACES it.  pr.jsonnet's tgm_fv_*_margin function DEFAULTS (x 2.5 / y 3 / z 5, zmin 3)
// are not the operating point: the production driver pdhd/wct-pr-perevt.jsonnet overrides
// them to x 2 / y 17.5 / z 18, i.e. the dvm margins (2.5 y, 3 z) PLUS 15 cm of flat shell,
// exactly as PDVD's driver does.  Carrying both the shell and this surface would count the
// 15 cm twice, so ONE knob moves both halves: curved_fv swaps the fiducial AND drops the
// y/z entries of fv_tolerance to the cushion alone (curved_fv_margin_y/z).  The consumer
// therefore passes cushion 0 here and keeps the cushion in fv_tolerance.
//   Consequence, measured (doc pdhd/09 sec 9.5, doc pdvd/43 sec 6.2, compared in doc
// pdvd/49 sec 5): PDHD's p90 surface sits INSIDE the shell it replaced over ~76 % of the
// drift, PDVD's over only ~25 % (y) / 44 % (z).  So the same knob moves the TOTAL TGM
// count in opposite directions -- PDHD 1478 -> 1664 (+12.6 %), PDVD 2148 -> 2095 (-2.5 %)
// -- while LONG-TRACK (> 2 m) TGM rises on BOTH, 430 -> 464 and 754 -> 769.  The opposite
// total signs are a short-cluster effect (TaggerCheckTGM.cxx:1066), not a disagreement
// about the boundary.  NB the often-quoted PDVD "-33 %" is the WITHDRAWN d50 arm
// (doc pdvd/41 sec 13.3), not production.
//
// The polygon spans BOTH drift volumes and is continuous across the cathode, exactly as
// pdhd_pr_fv's box is, so a cathode-crossing track is not an "exiter" at x = 0.  In x the
// cushion insets the two anode faces only, for the same reason.
//
// NOTHING IMPORTS THIS FILE -- that is its byte-identity argument, the same one
// protodunevd/curved_fiducial.jsonnet and crp_gap_fiducial.jsonnet make.  Wiring it into
// pr.jsonnet's tagger fiducial is a separate, knobbed change with an A/B of its own.

local wc = import 'wirecell.jsonnet';

// --- the volume pdhd_pr_fv spans, cfg/pgrapher/experiment/pdhd/pr.jsonnet:1258-1261 (cm) ---
local XW = 357.985;    // anode face, both drifts
local YLO = 7.61;      // low-y wall   (PDVD: -YW)
local YHI = 606.0;     // high-y wall  (PDVD: +YW)
local ZLO = 0.234345;  // upstream z wall
local ZHI = 462.297;   // downstream z wall
local CATH = 2.54;     // cathode face (per-face FV_xmax, clus.jsonnet:84)

// {dc: inset at the cathode face (cm), x1: end of the plateau, x2: foot of the ramp}.
// dc = 0 => that wall is at the nominal wall at every x.  All eight default to 0: PDHD
// has no measured trapezoid surface (see note 2 above); pass `profile` for a real one.
function(yp_g02={dc: 0.00, x1: CATH, x2: CATH},
         yp_g13={dc: 0.00, x1: CATH, x2: CATH},
         ym_g02={dc: 0.00, x1: CATH, x2: CATH},
         ym_g13={dc: 0.00, x1: CATH, x2: CATH},
         zm_g02={dc: 0.00, x1: CATH, x2: CATH},
         zm_g13={dc: 0.00, x1: CATH, x2: CATH},
         zp_g02={dc: 0.00, x1: CATH, x2: CATH},
         zp_g13={dc: 0.00, x1: CATH, x2: CATH},
         cushion_x=0.0, cushion_y=0.0, cushion_z=0.0,  // cm, inward from the surface
         name_prefix='pdhdcurved',
         // profile (doc pdhd/09): an object with the same eight keys whose values are
         // explicit knot lists [[|x|, inset], ...] in cm, anode face -> cathode face,
         // e.g. one entry of curved_fiducial_profiles.jsonnet (the exit-gap p80 / p90
         // surfaces).  null (the default) => the eight trapezoids above, i.e. flat.
         profile=null)

  local W = if profile == null
            then { yp_g02: yp_g02, yp_g13: yp_g13, ym_g02: ym_g02, ym_g13: ym_g13,
                   zm_g02: zm_g02, zm_g13: zm_g13, zp_g02: zp_g02, zp_g13: zp_g13 }
            else profile;

  // knots of one wall of one volume, anode face -> cathode face, as [|x|, inset] in cm:
  // either the trapezoid's own corners, or an explicit knot list passed through.
  local knots(p) =
    if std.isArray(p) then p else
    [[XW, 0.0]]
    + (if p.x2 < XW then [[p.x2, 0.0]] else [])
    + (if p.x1 > CATH then [[p.x1, p.dc]] else [])
    + [[CATH, p.dc]];

  // ... placed on a wall of one drift volume.  sgn = -1 for g02 (x<0, APA0/2), +1 for
  // g13 (x>0, APA1/3); inward = +1 for a low wall (interior at larger coordinate),
  // -1 for a high wall.
  local side(p, sgn, wallpos, inward, cush) =
    [[sgn * std.min(k[0], XW - cushion_x), wallpos + inward * (k[1] + cush)] for k in knots(p)];

  local xy =
    side(W.ym_g02, -1, YLO, 1, cushion_y)
    + std.reverse(side(W.ym_g13, 1, YLO, 1, cushion_y))
    + side(W.yp_g13, 1, YHI, -1, cushion_y)
    + std.reverse(side(W.yp_g02, -1, YHI, -1, cushion_y));

  local xz =
    side(W.zm_g02, -1, ZLO, 1, cushion_z)
    + std.reverse(side(W.zm_g13, 1, ZLO, 1, cushion_z))
    + side(W.zp_g13, 1, ZHI, -1, cushion_z)
    + std.reverse(side(W.zp_g02, -1, ZHI, -1, cushion_z));

  // drop a repeated vertex (a wall with x2 == XW has no ramp foot to name)
  local dedup(pts) = [pts[i] for i in std.range(0, std.length(pts) - 1)
                      if i == 0 || pts[i][0] != pts[i - 1][0] || pts[i][1] != pts[i - 1][1]];

  local XYP = dedup(xy);
  local XZP = dedup(xz);

  // PolyFiducial stacks polygonal slabs along `axis` and the corners are the two
  // TRANSVERSE coordinates in the order (axis+1, axis+2) mod 3 (PolyFiducial.cxx:62-65):
  //   axis 2 (z) -> corners are (x, y);  axis 1 (y) -> corners are (z, x).
  // One slab each, spanning the whole detector along its axis: the axis span is only a
  // bounding-box short circuit (PolyFiducial.cxx:139 returns false outside m_bb before
  // any slab is tried), the transverse polygon is the cut, and the two planes are
  // combined by the composite below -- slabs INSIDE one PolyFiducial are OR-ed, so the
  // AND has to be the composite's.
  local span = 1000.0;

  local xy_poly = {
    type: 'PolyFiducial',
    name: name_prefix + '-xy',
    data: {
      axis: 2,
      slabs: [{ min: -span * wc.cm, max: span * wc.cm,
                corners: [[p[0] * wc.cm, p[1] * wc.cm] for p in XYP] }],
    },
  };
  local xz_poly = {
    type: 'PolyFiducial',
    name: name_prefix + '-xz',
    data: {
      axis: 1,
      slabs: [{ min: -span * wc.cm, max: span * wc.cm,
                corners: [[p[1] * wc.cm, p[0] * wc.cm] for p in XZP] }],
    },
  };
  local composite = {
    type: 'CompositeFiducial',
    name: name_prefix + '-fv',
    data: {
      logic: 'and',
      fiducials: [wc.tn(xy_poly), wc.tn(xz_poly)],
    },
  };

  {
    polys: [xy_poly, xz_poly],
    composite: composite,
    tn: wc.tn(composite),       // reference the volume by this type:name
    configs: [xy_poly, xz_poly, composite],
    // the vertex lists in cm, for a doc or a plot
    boundary_xy: XYP,
    boundary_xz: XZP,
  }
