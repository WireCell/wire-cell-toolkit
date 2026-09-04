// PDVD CRP / CRU / cathode structure-exclusion fiducial volume.
//
// Counterpart of cfg/pgrapher/experiment/sbnd/cathode_fiducial.jsonnet.  SBND has
// exactly ONE structural gap (the CPA) and models it from engineering drawings;
// PDVD has four families of gap and this file models them from MEASUREMENT, for
// the reason doc pdvd/35 sec 3 sets out: the width that matters to a fiducial
// volume is the width the reconstruction loses, and on PDVD that is not the
// width the geometry declares -- in one direction for the cathode and in the
// other for the CRU seams.
//
// PDVD geometry, read off the 16 AnodePlane 'sensvol' boxes (doc 35 sec 2):
//   2 drift volumes; 2 CRPs per drift (y<0 and y>0), each 3357.8 x 2976.2 mm
//   ~ the 3.0 x 3.4 m^2 ProtoDUNE-VD CRP; each CRP is 4 CRUs, 2 in y x 2 in z.
// So the seams are three different things and must not be given one width:
//   y = 0        CRP <-> CRP   (anodes 0,1 <-> 2,3)   12.2 mm mechanical
//   |y| = 168.5  CRU <-> CRU   (face 0 <-> face 1)     0.2 mm mechanical
//   z = 149.65   CRU <-> CRU   (anode 0 <-> anode 1)   2.6 mm bot / 1.0 mm top
//   |x| < 3.0    the cathode slab                     60.0 mm mechanical
//
// Measured half-widths (doc 35 sec 3; fv_gap_measure.py A2 = the median empty
// interval a crossing track loses, minus the same-axis control's point-pitch
// floor, over the 120-event d28dlfp arm):
//   cathode      4.08 cm   vs 3.00 geometric  -> the ONLY family that dilates
//   y = 0        0.61 cm   vs 0.61 geometric  -> exactly its mechanical gap
//   |y| = 168.5  0.03-0.15 vs 0.01 geometric  -> at the measurement floor
//   z = 149.65   0.00-0.06 vs 0.13/0.05 geom  -> consistent with zero
// The CRU boxes are therefore built at their GEOMETRIC width: nothing was
// measured that would justify dilating them, and doc 35 sec 3 says so.
//
// Returns BoxFiducial + CompositeFiducial{logic:'or'} component configs -- the
// same generic IFiducial primitives SBND uses, so this is a point-in-structure
// test reusable anywhere.  Returns:
//   { boxes:     [<8 BoxFiducial cfgs>],
//     composite: <CompositeFiducial cfg, OR over all boxes>,
//     tn:        'CompositeFiducial:<prefix>-gaps',   // reference by this tn
//     configs:   boxes + [composite] }                // inject this whole list
//
// NOTE (doc 35 sec 6): nothing imports this file yet.  QLMatching's
// `cathode_fiducial` is read at QLMatching.cxx:5149 for the cathode-end
// at_x_boundary flag only, where PDVD's flat 6 cm slab gains nothing from a 3-D
// test and the CRU seams are invisible.  The consumer this was built for is a
// FiducialUtils `structure_fiducial` -- a step inside it scored TRANSPARENT
// rather than dead in check_dead_volume / check_signal_processing -- which is
// its own round.

local wc = import 'wirecell.jsonnet';

// Cushions default to ZERO, unlike SBND's 0.5 cm, and the difference is not an
// oversight.  SBND dilates a MECHANICAL model (pad / tube / knuckle reaches off
// a drawing) to cover whatever the reconstruction adds on top.  The two widths
// that matter here were measured on the reconstruction's own output, so they
// already contain that; adding a cushion would double-count it.  The remaining
// families are at geometric width precisely because no excess was measurable,
// which is not a licence to pad them either.  A5 of doc 35 backs this from the
// other side: at the outer y/z walls the imaged-point density goes from 0 to
// full within ~2 mm, so there is no soft edge for a cushion to cover.  The
// arguments stay exposed for a consumer that wants margin of its own.
function(cx=0.0 * wc.cm, cy=0.0 * wc.cm, cz=0.0 * wc.cm,
         cathode_reach=4.075 * wc.cm,      // MEASURED (geometric 3.0)
         crp_y0_reach=0.610 * wc.cm,       // measured 0.61 == geometric 0.61
         cru_y_reach=0.010 * wc.cm,        // geometric; measured 0.03-0.15, at the floor
         cru_z_reach_bot=0.130 * wc.cm,    // geometric; measured ~0.06
         cru_z_reach_top=0.050 * wc.cm,    // geometric; measured ~0.00
         include_cru=true,                 // false => cathode + the CRP seam only
         name_prefix='pdvdgap')

  // --- envelope, from the sensvol union (cm) ---
  local XW = 339.91 * wc.cm;    // shield plane, both drifts
  local YW = 336.39 * wc.cm;    // outer y wall
  local ZLO = 0.813 * wc.cm;    // upstream z wall
  local ZHI = 298.435 * wc.cm;  // downstream z wall
  local CRU_Y = 168.50 * wc.cm; // the two CRU seams in y, +- this
  local CRU_Z = 149.65 * wc.cm; // the CRU seam in z (same centre both drifts)

  // full-extent spans, dilated by the same cushion so the slabs reach the walls
  local xfull = { lo: -(XW + cx), hi: XW + cx };
  local yfull = { lo: -(YW + cy), hi: YW + cy };
  local zfull = { lo: ZLO - cz, hi: ZHI + cz };

  local slab(name, x, y, z) = { name: name, x: x, y: y, z: z };
  local band(centre, reach, cush) = { lo: centre - (reach + cush), hi: centre + (reach + cush) };

  local cathode = [
    slab('cathode', band(0, cathode_reach, cx), yfull, zfull),
  ];
  local crp_seam = [
    slab('crp_y0', xfull, band(0, crp_y0_reach, cy), zfull),
  ];
  local cru_seams = [
    slab('cru_y_pos', xfull, band(CRU_Y, cru_y_reach, cy), zfull),
    slab('cru_y_neg', xfull, band(-CRU_Y, cru_y_reach, cy), zfull),
    // the z seam is 2.6 mm in the bottom drift and 1.0 mm in the top, so it is
    // two boxes, split at the cathode, not one.
    slab('cru_z_bot', { lo: xfull.lo, hi: 0.0 }, yfull, band(CRU_Z, cru_z_reach_bot, cz)),
    slab('cru_z_top', { lo: 0.0, hi: xfull.hi }, yfull, band(CRU_Z, cru_z_reach_top, cz)),
  ];

  local slabs = cathode + crp_seam + (if include_cru then cru_seams else []);

  local boxes = [
    { type: 'BoxFiducial',
      name: '%s-%s' % [name_prefix, s.name],
      data: { bounds: { tail: { x: s.x.lo, y: s.y.lo, z: s.z.lo },
                        head: { x: s.x.hi, y: s.y.hi, z: s.z.hi } } } }
    for s in slabs ];

  local composite = {
    type: 'CompositeFiducial',
    name: name_prefix + '-gaps',
    data: { logic: 'or', fiducials: ['BoxFiducial:' + b.name for b in boxes] },
  };

  {
    boxes: boxes,
    composite: composite,
    tn: 'CompositeFiducial:' + name_prefix + '-gaps',
    configs: boxes + [composite],
  }
