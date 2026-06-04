// Hand-declared dead regions for SBND (beyond the chndb "bad" channel list).
//
// Event-1720 region: the W collection wires are dead and the U/V signal is
// locally distorted near y~0, z~251 cm, leaving a ~4 cm gap that fragments a
// single through-going cosmic inside examine_bundles' relaxed connectivity
// graph.  We declare that spot dead in two complementary ways, driven from ONE
// channel list so the two can never drift apart:
//   - dead WINDS (clus/PointTreeBuilding "inject_dead_winds"): the FUNCTIONAL fix
//     -- makes Grouping::is_good_point()/test_good_point() (hence the relaxed-graph
//     bridge in connect_graph_relaxed) treat the gap as a crossable dead region.
//   - dead BLOB  (img/MaskSlice "masked_channels"): VISUAL only -- so the region
//     shows up in the Bee deadarea display.  Inert for the split.
//
// Channels are half-open [lo, hi): W [4800,4806) = exactly the 6 chndb-bad W
// channels (4800-4805).  TPC1 (anode ident 1) mirrors TPC0 by the +5638 per-APA
// channel offset (W [10438,10444) is already chndb-bad -- confirms the mirror).
{
    local apa_nchan = 5638,

    // TPC0 (anode ident 0) channel idents of the bad region.
    local tpc0_channels =
        std.range(988, 990)        // U [988, 991)
        + std.range(2970, 2983)    // V [2970, 2984)
        + std.range(4800, 4805),   // W [4800, 4806)  (already chndb-bad)

    // Per-anode channel list (TPC0 set shifted by anode_ident * 5638).
    // Used by the imaging masked_channels dead-blob branch.
    channels(anode_ident):: [c + anode_ident * apa_nchan for c in tpc0_channels],

    // Per-anode dead-wind region: channels + full-drift x window.  xbeg/xend are
    // in WCT system-of-units length; +-1e9 comfortably exceeds any SBND drift
    // coordinate (|x| < ~2.1 m), i.e. "all drift" (the defect is a fixed
    // anode-plane region, not time-localized).  Used by inject_dead_winds.
    region(anode_ident):: {
        channels: $.channels(anode_ident),
        xbeg: -1e9,
        xend: 1e9,
    },
}
