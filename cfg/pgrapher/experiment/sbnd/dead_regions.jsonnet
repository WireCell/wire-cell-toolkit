// Hand-declared dead regions for SBND (beyond the chndb "bad" channel list).
//
// The SBND anode has a fixed defect band -- the "middle 6-ch W-plane dead area"
// (TPC0 W channels 4800-4805): the collection (W) plane has no wires and the
// induction (U/V) signal is locally distorted.  The gap is UNIFORM top-to-bottom
// (a full vertical column where W is dead), so it must be treated as a dead
// region for all three planes over the WHOLE column, not just the y~0 center.
//
// We declare the W channels dead and flag the region as a "gap": the functional
// fix is in clus (PointTreeBuilding "inject_dead_winds" -> Grouping dead-gap
// registry).  Because a W wind maps to a single z (collection wires are vertical),
// keying the gap off the W channels makes the whole vertical column crossable in
// Grouping::is_good_point()/test_good_point() (hence the relaxed-graph bridge in
// connect_graph_relaxed) WITHOUT injecting any U/V wires -- a U/V wire is diagonal
// in (y,z), so declaring it dead would pollute the entire diagonal far from the gap.
//
// Channels are half-open [lo, hi): W [4800,4806) = exactly the 6 chndb-bad W
// channels (4800-4805).  TPC1 (anode ident 1) mirrors TPC0 by the +5638 per-APA
// channel offset (W [10438,10444) is already chndb-bad -- confirms the mirror).
{
    local apa_nchan = 5638,

    // TPC0 (anode ident 0) W channel idents of the bad region.
    local tpc0_w_channels = std.range(4800, 4805),   // W [4800, 4806)  (already chndb-bad)

    // Per-anode dead-gap region: the W channels + full-drift x window, flagged as
    // a "gap".  xbeg/xend are in WCT system-of-units length; +-1e9 comfortably
    // exceeds any SBND drift coordinate (|x| < ~2.1 m), i.e. "all drift" (the
    // defect is a fixed anode-plane region, not time-localized).  Used by
    // inject_dead_winds; gap:true routes the W winds into the dead-gap registry.
    region(anode_ident):: {
        channels: [c + anode_ident * apa_nchan for c in tpc0_w_channels],
        xbeg: -1e9,
        xend: 1e9,
        gap: true,
    },
}
