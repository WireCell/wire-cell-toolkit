// PDVD space-charge (curved) fiducial volume, MEASURED.
//
// The counterpart of the MicroBooNE prototype's ToyFiducial polygons
// (prototype_base/wire-cell/pid/src/ToyFiducial.cxx:60-135, 204-215): in each of the
// X-Y and X-Z planes the apparent detector boundary is not a box but a surface that is
// at the nominal wall near the anode and moves inward toward the cathode, and the two
// planes are tested INDEPENDENTLY and AND-ed (uBooNE: pnpoly(XY) && pnpoly(XZ), here
// CompositeFiducial{logic:'and'} of two PolyFiducial).
//
// PDVD's numbers below are measured on data, not adopted: doc pdvd/41 maps the apparent
// wall position against drift x from 120 Q/L-matched cosmic events (7.8 M imaged points
// with a cluster T0), per wall and per drift volume, with an anode-side control that
// returns the geometric wall to within 1.1 cm and a straight-line extrapolation test
// showing the effect is charge DISPLACEMENT, not charge loss.  See doc 41 sec 9 for the
// fit table and sec 5-7 for the map.
//
// Each wall of each drift volume is one trapezoid in |x| (cm):
//
//     inset(|x|) = dc                          for |x| <= x1     (the cathode plateau)
//                = dc * (x2 - |x|)/(x2 - x1)   for x1 < |x| < x2 (the ramp)
//                = 0                           for |x| >= x2     (at the nominal wall)
//
// x1 = the cathode face (3 cm) recovers uBooNE's two-segment shape exactly; x2 <= the
// anode face is REQUIRED, because the displacement a drifting charge accumulates is
// zero at the anode by construction.
//
// The eight profiles are NOT related by symmetry: doc 41 sec 6.1 rejects the bottom/top
// mirror (chi2 151/34 on z+) and the y+/y- mirror (110/34) on this data set, while the
// Y-vs-X and Z-vs-X factorization -- the thing that makes two polygons enough -- holds
// (chi2/ndf <= 1.3 in z slabs and CRU quarters).  So all eight are carried separately
// and there is one slab per plane.
//
// CUSHION.  These vertices are the calibrated surface itself, cushion 0.  MicroBooNE
// applies its cushion two different ways and this file supports both:
//   - the containment / physics fiducial insets the polygon VERTICES by
//     boundary_dis_cut = 3 cm (4 cm on the z walls: 3 + a hard-coded 1) at every wall
//     and both x faces -- pass cushion_y / cushion_z / cushion_x here;
//   - the cosmic tagger keeps the uncushioned surface and probes it with a tolerance
//     BAND of 2.0-2.8 cm (Cosmic_tagger.h:34-36 stm/tgm/tgm-no-flash vectors) -- in WCT
//     that is the taggers' fv_tolerance, which offsets the POINT instead of the vertex
//     (FiducialUtils.cxx:79-119).  The two agree to cos(ramp angle) = 0.989 here (the
//     steepest ramp is y- bottom, 9.2 cm over 61 cm), i.e. under 4 mm on a 3 cm cushion.
// PDVD's tagger margins today are 2.5 cm (x, y) and 3 cm (z) plus a flat 15 cm
// space-charge allowance in y and z (doc pdvd/35); this surface REPLACES the 15 cm, so
// a consumer wanting uBooNE's arrangement passes cushion 0 here and drops fv_tolerance
// back to the 2.5 / 2.5 / 3 cm cushion.
//
// The polygon spans BOTH drift volumes and is continuous across the 6 cm cathode slab,
// exactly as pr.jsonnet's pdvd_pr_fv box is, so a cathode-crossing track is not an
// "exiter" at x = 0.  In x the cushion insets the two anode faces only, for the same
// reason.
//
// NOTHING IMPORTS THIS FILE -- that is its byte-identity argument, the same one
// crp_gap_fiducial.jsonnet makes.  Wiring it into pr.jsonnet's tagger fiducial is a
// separate, knobbed change with an A/B of its own (doc 41 sec 9.5).

local wc = import 'wirecell.jsonnet';

// --- the sensitive-volume envelope, from the 16 AnodePlane 'sensvol' boxes (cm) ---
local XW = 339.91;    // anode (shield) plane, both drifts
local YW = 336.39;    // outer y wall
local ZLO = 0.813;    // upstream z wall
local ZHI = 298.435;  // downstream z wall
local CATH = 3.0;     // cathode face

// doc 41 sec 9.1, d50 half-density estimator, 20 cm drift bins, 200 event bootstraps:
// {dc: inset at the cathode face (cm), x1: end of the plateau, x2: foot of the ramp}.
// dc = 0 => that wall is at the nominal wall at every x (y+ bottom: 0.15 +- 0.69 cm).
function(yp_bot={dc: 0.00, x1: CATH, x2: CATH},        // flat: amplitude under 2 sigma
         yp_top={dc: 2.76, x1: CATH, x2: 126.82},
         ym_bot={dc: 9.22, x1: 114.79, x2: 176.19},    // the one profile a plateau fits better
         ym_top={dc: 2.45, x1: CATH, x2: XW},
         zm_bot={dc: 11.15, x1: CATH, x2: 200.95},
         zm_top={dc: 3.78, x1: CATH, x2: XW},
         zp_bot={dc: 17.66, x1: CATH, x2: 205.14},
         zp_top={dc: 11.02, x1: CATH, x2: 164.61},
         cushion_x=0.0, cushion_y=0.0, cushion_z=0.0,  // cm, inward from the surface
         name_prefix='pdvdcurved')

  // knots of one wall of one volume, anode face -> cathode face, as [|x|, inset] in cm.
  local knots(p) =
    [[XW, 0.0]]
    + (if p.x2 < XW then [[p.x2, 0.0]] else [])
    + (if p.x1 > CATH then [[p.x1, p.dc]] else [])
    + [[CATH, p.dc]];

  // ... placed on a wall of one drift volume.  sgn = -1 bottom (x<0), +1 top (x>0);
  // inward = +1 for a low wall (interior at larger coordinate), -1 for a high wall.
  local side(p, sgn, wallpos, inward, cush) =
    [[sgn * std.min(k[0], XW - cushion_x), wallpos + inward * (k[1] + cush)] for k in knots(p)];

  local xy =
    side(ym_bot, -1, -YW, 1, cushion_y)
    + std.reverse(side(ym_top, 1, -YW, 1, cushion_y))
    + side(yp_top, 1, YW, -1, cushion_y)
    + std.reverse(side(yp_bot, -1, YW, -1, cushion_y));

  local xz =
    side(zm_bot, -1, ZLO, 1, cushion_z)
    + std.reverse(side(zm_top, 1, ZLO, 1, cushion_z))
    + side(zp_top, 1, ZHI, -1, cushion_z)
    + std.reverse(side(zp_bot, -1, ZHI, -1, cushion_z));

  // drop a repeated vertex (a wall with x2 == XW has no ramp foot to name)
  local dedup(pts) = [pts[i] for i in std.range(0, std.length(pts) - 1)
                      if i == 0 || pts[i][0] != pts[i - 1][0] || pts[i][1] != pts[i - 1][1]];

  local XYP = dedup(xy);
  local XZP = dedup(xz);

  // PolyFiducial stacks polygonal slabs along `axis` and the corners are the two
  // TRANSVERSE coordinates in the order (axis+1, axis+2) mod 3 (PolyFiducial.cxx:62-65):
  //   axis 2 (z) -> corners are (x, y);  axis 1 (y) -> corners are (z, x).
  // One slab each, spanning the whole detector along its axis: the axis span is only a
  // bounding-box short circuit here, the transverse polygon is the cut, and the two
  // planes are combined by the composite below.
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
