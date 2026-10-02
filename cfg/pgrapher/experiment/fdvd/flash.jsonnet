// FD-VD OpFlashFinder data for the in-toolkit low-energy reconstruction (fdvd_sim doc 16).
//
// Exactly the FROZEN doc 04 point the python Q-L chain of fdvd_sim docs 05-15 consumed:
// fdvd_sim/output/fdvd04_flash/lateoff_tw10_thr2.5/config.json (wcp WireCell/wcp-porting-validation),
// written by the doc 04 sweep (stageB/sweep_flash.py).  It differs from fdvd_sim/cfg/fdvd-flash.jsonnet
// in min_fired_pds 0 (there 1) and min_total_pe 0 (there 10): the sweep left the flash list uncut and
// applied its cuts in scoring, and the Q-L matcher's own group cut (>= 3 OpDets with >= 1 PE,
// ql_purity.PMIN_STORE) supersedes the 1-PD cut.  Keeping 0 / 0 is what makes the in-toolkit flashes
// hash-identical to the python reference (doc 16 gate H).
//   remove_late_light false   every late-light setting loses MARLEY flashes in mixed events (doc 04 sec 1)
//   flash_threshold 2.5, tail merge on (window 10 us)
// FD-VD geometry: 184 X-ARAPUCAs (fdvd-opdet-geom.json), 368 OpChannels -> 184 (fdvd-opch-map.json).
{
  nchan: 184,
  geom_file: 'pgrapher/experiment/fdvd/fdvd-opdet-geom.json',
  channel_map_file: 'pgrapher/experiment/fdvd/fdvd-opch-map.json',
  group_by_side: false,     // must stay off: it splits at x = 0 (cathode PDs vs walls)
  flash_refine: false,      // must stay off: its adjacency grid is the PDHD layout
  offset_us: 0,
  remove_late_light: false,
  flash_threshold: 2.5,
  flash_tail_merge: true,
  tail_window_us: 10,
  tail_min_width_us: 1,
  tail_pe_frac: 0.7,
  tail_pe_ratio: 1,
  min_fired_pds: 0,
  min_total_pe: 0,
  min_fired_pe: 1,
}
