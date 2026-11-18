# Truth-match id conventions (SBND `truth_match`)

With `truth_match` on (MC only), `SbndPrMagnifyTrackingVisitor` writes the truth-matching trees to `tracking-pr.root`: `T_truth_track`, `T_overlap`, `T_match` and `T_pf`, plus `T_blob_truth` when `blob_truth` is on. The code is in `root/src/SbndPrTruthMatch.cxx`.

The truth ids in these trees come from the step-1 labeler, larwirecell `wclsTensorSetLabeler` with `truth_detail` (introduced in larwirecell [`73c79bc`](https://github.com/HaiwangYu/larwirecell/commit/73c79bc48de7bd631f1d096b96512a41aa161de7)). It writes each blob's dominant `trackid` and the `truth_track` tensor into the step-1 output, and step 2 reads them. This page gives the id conventions. The full branch tables are in the [issue-41 reference](https://github.com/HaiwangYu/wire-cell-toolkit-ai-helper/blob/0898fb8/issues/41-sbnd-reco-truth-matching/41-sbnd-reco-truth-matching.md#8-reference-the-new-trees).

## Blob `trackid`: the dominant track

1. Each SimEnergyDeposit is projected to wires and ticks after SCE, drift with the deposit's own time, diffusion and the SP smearing.
2. A blob accepts a deposit within its wire and tick bounds plus slop (1 wire, 2 ticks) and 3σ.
3. The accepted deposits are weighted by `NumElectrons`, or by `Energy()` if that is 0, and summed per track.
4. The largest sum is the blob's `trackid`, and the runner-up is `trackid2`.
5. `truth_q`, `truth_q_dom` and `truth_q2` are the total, dominant and runner-up sums, in electrons.

## Sentinels

- **`trackid = -1`** is the only trackid sentinel: the blob accepted no deposit (a ghost or noise blob), and every blob on data.
  - Step 2 also maps a `trackid` with no `truth_track` row to -1, as a guard. With `truth_detail`, every dominant track has a row.
  - The points of such blobs are the "unlabeled" ones: `T_overlap.trackid = -1`, `T_match.n_unlabeled`.
- **`T_overlap.nu_idx`:** -1 = a labeled track not from a beam neutrino, -2 = an unlabeled row (`trackid = -1`).
- **`T_blob_truth.trackid2`:** -1, with `truth_q2` 0, when the blob has only one track.

## larg4 negative ids, and unsaved particles

- larg4 does not save every particle. A saved particle's deposits carry its own `TrackID`. Deposits of a descendant it did not save carry the **negative** of the saved ancestor's id. Examples are delta rays and EM shower daughters that were not kept.
- The labeler uses `abs(TrackID)`, so that charge counts under the saved ancestor everywhere: blob `trackid`, `q_readout`, `q_in_blobs`. No negative id reaches the trees.
- There is no special case by origin. Unsaved EM descended from a neutron (capture γ's and their electrons, for example) folds into its nearest saved ancestor, as any other unsaved descendant does.
- A saved particle keeps its own row, and its `nu_idx` is the interaction of its MCTruth association.

## `T_truth_track` grouping ids

One row per G4 track with charge inside the readout window, or dominating a blob.

| id | rule |
|---|---|
| `nu_idx` | 0-based beam-neutrino MCTruth key = `T_truth_nu.nu_idx`; -1 = not from a beam neutrino (cosmics and the like). The Bee `sed` sets use 1-based (0 = not a neutrino); the trees do not |
| `bee_trackid` | `trackid`, except that a Michel electron (e±, `Decay`, mother μ±) maps to its mother muon: the Bee truth-layer merge |
| `ancestor_trackid` | follow `Mother()` while the mother is a stored MCParticle; the generator particle or the cosmic primary |
| `prong_trackid` | the truth particle a reco particle is compared with. Starting from the track, move to the mother while the track is (a) an e±/γ whose mother is e±/γ (a shower is one prong), (b) an ionization electron (`eIoni`, `muIoni`, `hIoni`, `ionIoni`), or (c) a muon's radiative e±/γ (`muPairProd`, `muBrems`). A Michel electron stops the walk, so it stays its own prong |
| `process` | the CellTree G4 process code (larwirecell `aiml/G4ProcessCode.h`), with the synthetic Michel = 10001 and -1 for an unknown process |

## Reco side of a shower

- A segment can sit in two `PRShower`s; the first shower that claims it is the one used.
- When the segment's `T_kine` row is a shower row, its `particle_id` (`T_overlap`) and `shower_id` (`T_segment`) are that row's `kine_particle_id`, so they agree with `kine_index`.
