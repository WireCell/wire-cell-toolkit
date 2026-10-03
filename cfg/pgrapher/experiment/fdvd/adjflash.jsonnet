// FD-VD FdvdAdjOpFlashFinder data (AdjOpHits flash clustering) for lowe-reco.jsonnet flash_finder='adjophits'
// (fdvd_sim doc 21).
//
// The point is the duneana SolarNuAna in-module AdjOpHits setting (fdvd_sim doc 19 'sna_inmodule', doc 20 'adj0'):
// backward window 1.0 us, forward window 1.6 us, radius 600 cm, no NHit / PE / trigger-PE cuts, hit duplicates on,
// flash time = StartTime of the brightest member.  Doc 21 re-selected it on the tuning split of the doc 18 overlays
// (12 points of window x radius) by the rule pre-written in doc 21 sec 4.
// The C++ defaults are the official FD-VD SolarOpFlash setting (+-30 ns, 500 cm, NHit 3, PE 1.5); every key that
// differs is given here.
{
  nchan: 184,
  geom_file: 'pgrapher/experiment/fdvd/fdvd-opdet-geom.json',
  channel_map_file: 'pgrapher/experiment/fdvd/fdvd-opch-map.json',
  time_var: 'StartTime',
  min_time_us: 1.0,
  max_time_us: 1.6,
  radius_cm: 600.0,
  nhit: 0,
  pe: 0.0,
  trigger_pe: 0.0,
  hot_threshold: 0.3,
  hit_duplicates: true,
}
