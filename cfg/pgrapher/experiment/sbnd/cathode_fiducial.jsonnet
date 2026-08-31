// SBND cathode-plane (CPA) structure-exclusion fiducial volume.
//
// Faithful jsonnet port of sbnd_xin/sbnd_geometry/cathode_fiducial.py (the design
// + validation source of truth). The CPA is NOT a flat sheet: thin foil pads reach
// ~0.6 cm into the drift, the steel tube lattice ~2.7 cm, and the knuckle joints
// ~4.1 cm. We approximate the occupied region as a union of axis-aligned boxes,
// per TPC, in the wire-cell-toolkit frame (X=drift, cathode at X=0; TPC0 East X<0,
// TPC1 West X>0; Y vertical centered 0; Z beam, Z_toolkit = 250.5cm + Z_local).
//
// The whole region is dilated by an independent per-side cushion in each dimension
// (default 0.5 cm), applied to pad slab, tube bars and knuckles alike.
//
// Returns BoxFiducial + CompositeFiducial{logic:'or'} component configs (the toolkit
// IFiducial primitives), so this is a generic point-in-volume tool reusable anywhere
// (matching, pattern recognition, ...). Returns:
//   { boxes:     [<26 BoxFiducial cfgs>],
//     composite: <CompositeFiducial cfg, OR over all boxes>,
//     tn:        'CompositeFiducial:<prefix>-exclusion',   // reference by this tn
//     configs:   boxes + [composite] }                     // inject this whole list

local wc = import 'wirecell.jsonnet';

function(cx=0.5*wc.cm, cy=0.5*wc.cm, cz=0.5*wc.cm, pad=true, tube_hw=null,
         name_prefix='cpa')

  // --- structure constants (cm; from cathode_fiducial.py, mm->cm) ---
  local PAD_REACH  = 0.6*wc.cm;    // thin pad slab depth
  local TUBE_REACH = 2.7*wc.cm;    // pipe reach into drift
  local TUBE_HW    = 2.7*wc.cm;    // tube bar transverse half-width (pipe radius)
  local KNK_REACH  = 4.1*wc.cm;    // knuckle joint reach into drift
  local KNK_HY     = 5.0*wc.cm;    // knuckle half-size in y (100 mm tall)
  local KNK_HZ     = 6.5*wc.cm;    // knuckle half-size in z (130 mm)
  local PLANE_YH   = 207.0*wc.cm;  // bar long-extent (half) in y
  local PLANE_ZH   = 252.8*wc.cm;  // bar long-extent (half) in z
  local Z_CENTER   = 250.5*wc.cm;  // volTPCActive local z=0 -> toolkit z

  local hw = if tube_hw == null then TUBE_HW else tube_hw;

  // local z -> toolkit z
  local zc(z_local) = Z_CENTER + z_local;

  // X span (per TPC) for a feature reaching `reach` into the drift, + cushion cx.
  // TPC0 -> [-(reach+cx), 0]; TPC1 -> [0, +(reach+cx)]. (PER_TPC_REACH in the python
  // is empty -- both TPCs mirror-symmetric; this is the per-TPC asymmetry hook.)
  local xspan(tpc, reach) =
    local d = reach + cx;
    if tpc == 0 then { lo: -d, hi: 0.0 } else { lo: 0.0, hi: d };

  // lattice lines (label matches the python %+05.0f mm name; value in toolkit/local cm)
  local HORIZ_Y = [  // tubes running along z, at y gap-lines
    { lab: '+0000', v: 0.0*wc.cm }, { lab: '+1035', v: 103.5*wc.cm },
    { lab: '-1035', v: -103.5*wc.cm }, { lab: '+2070', v: 207.0*wc.cm },
    { lab: '-2070', v: -207.0*wc.cm } ];
  local VERT_Z = [   // tubes running along y, at z gap-lines (local z)
    { lab: '+0000', v: 0.0*wc.cm }, { lab: '+1275', v: 127.5*wc.cm },
    { lab: '-1275', v: -127.5*wc.cm } ];
  local KNK_Y = [    // knuckles, all on z=0 (local)
    { lab: '+0518', v: 51.75*wc.cm }, { lab: '-0518', v: -51.75*wc.cm },
    { lab: '+1552', v: 155.25*wc.cm }, { lab: '-1552', v: -155.25*wc.cm } ];

  // per-TPC box list (deep boxes first, then the bbox-derived pad slab)
  local boxes_for(tpc) =
    local xs_tube = xspan(tpc, TUBE_REACH);
    local xs_knk  = xspan(tpc, KNK_REACH);
    local yplane = PLANE_YH + cy;
    local zlo = zc(-PLANE_ZH) - cz;
    local zhi = zc(PLANE_ZH) + cz;
    local htubes = [
      { name: 'htube_y' + e.lab,
        tail: { x: xs_tube.lo, y: e.v - (hw + cy), z: zlo },
        head: { x: xs_tube.hi, y: e.v + (hw + cy), z: zhi } }
      for e in HORIZ_Y ];
    local vtubes = [
      { name: 'vtube_z' + e.lab,
        tail: { x: xs_tube.lo, y: -yplane, z: zc(e.v) - (hw + cz) },
        head: { x: xs_tube.hi, y: yplane,  z: zc(e.v) + (hw + cz) } }
      for e in VERT_Z ];
    local knuckles = [
      { name: 'knuckle_y' + e.lab,
        tail: { x: xs_knk.lo, y: e.v - (KNK_HY + cy), z: zc(-KNK_HZ) - cz },
        head: { x: xs_knk.hi, y: e.v + (KNK_HY + cy), z: zc(KNK_HZ) + cz } }
      for e in KNK_Y ];
    local deep = htubes + vtubes + knuckles;
    // thin full-plane pad slab: its Y/Z extent is the bounding box of the whole
    // tube/knuckle lattice (cushions included), so it covers the cathode plane out
    // to and including the edge tube bars.
    local pad_box =
      local xs_pad = xspan(tpc, PAD_REACH);
      local ymin = std.foldl(function(a, b) std.min(a, b.tail.y), deep, deep[0].tail.y);
      local ymax = std.foldl(function(a, b) std.max(a, b.head.y), deep, deep[0].head.y);
      local zmin = std.foldl(function(a, b) std.min(a, b.tail.z), deep, deep[0].tail.z);
      local zmax = std.foldl(function(a, b) std.max(a, b.head.z), deep, deep[0].head.z);
      [{ name: 'pad', tail: { x: xs_pad.lo, y: ymin, z: zmin },
                      head: { x: xs_pad.hi, y: ymax, z: zmax } }];
    (if pad then pad_box else []) + deep;

  local all = [{ tpc: t, b: bx } for t in [0, 1] for bx in boxes_for(t)];
  local boxes = [
    { type: 'BoxFiducial',
      name: '%s-tpc%d-%s' % [name_prefix, e.tpc, e.b.name],
      data: { bounds: { tail: e.b.tail, head: e.b.head } } }
    for e in all ];
  local composite = {
    type: 'CompositeFiducial',
    name: name_prefix + '-exclusion',
    data: { logic: 'or', fiducials: ['BoxFiducial:' + b.name for b in boxes] },
  };

  {
    boxes: boxes,
    composite: composite,
    tn: 'CompositeFiducial:' + name_prefix + '-exclusion',
    configs: boxes + [composite],
  }
