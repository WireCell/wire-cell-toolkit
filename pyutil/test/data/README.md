# pyutil test data

- `dv2max_golden.json` (about 0.5 MB): 26 held-out DV2-MAX calls (fold 0 of the r7 DL-vertex dump, which the
  `r7_cloud` model was not trained on) and 3 synthetic calls with fewer than 5 points, with the answer the code
  gave when the file was made: `[x, y, z, P(<1 cm)]` in the raw reco frame, or empty for no prediction.
  `../test_dv2max.py` replays all calls and `../doctest_dv2_vertex.cxx` replays the first ones through the C++
  wrapper; both need the bundle `dv2max/sbnd/r7_cloud/bundle.json` (wire-cell-data) and skip without it.
  Calls are chosen so that the top-1 candidate wins by a clear ranker-score margin, so the answers do not depend
  on floating-point details of the machine.  Made by `../make_dv2max_golden.py`; `meta` records the inputs, the
  bundle sha256, the library versions and the tolerances (0.01 cm, 0.001 in score).  Regenerate it when the
  model bundle or the inference code changes on purpose.
- `dv2max_space_golden.json` (about 0.35 MB): 10 held-out calls of the space-point bundle
  `dv2max/sbnd/aurora_space/bundle.json` (its hold-out, fold 0 of the Aurora `aurora-r2` production, which the
  model was not trained on): BNB and nue, both passes, two without a traditional vertex, each with its cloud and
  space points (base64 of little-endian float32; cluster ids int32) and the answer the code gave when the file was
  made. `../test_dv2max_space.py` (`TestSpaceGolden`) replays them and skips without the bundle. Same selection
  idea (top-1 wins by a ranker-score margin of at least 0.5); tolerances 0.001 cm and 0.001 in score. Made by
  `../make_dv2max_space_golden.py`; the header records the inputs, the bundle sha256 and the library versions.
