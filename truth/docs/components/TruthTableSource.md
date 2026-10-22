# TruthTableSource

Serves the particle tables of one or more truth files (`truth/python/dump_truth_table.py`, written inside Stage A) as tensor sets: one set per configured ident, then EOS. The two tensors have the layout the SBND `wclsTruthInformationAttacher` produced, so `MultiAlgBlobClustering` (`aux_datatypes`) and `SbndPrMagnifyTrackingVisitor` consume them unchanged.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TruthTableSource` |
| Concrete class | `WireCell::Truth::TruthTableSource` |
| Node category | source |
| Primary interface | `ITensorSetSource` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inname` | truth file name; `%d` is replaced by the ident (default `truth_%d.npz`) |
| `idents` | the idents to serve, in order (default `[0]`) |
| `prefix` | datapath prefix of the tensors, `%d` = ident (default `truth/%d`) |
| `stop_on_missing` | a missing file ends the stream instead of throwing (default false) |

## Output tensors

| datatype | datapath | shape | columns (`md["columns"]`) |
| --- | --- | --- | --- |
| `truth_nu` | `<prefix>/nu` | double (Nn, 12) | `nu_idx, pdg, ccnc, mode, int_type, flavor, E, vtx_x, vtx_y, vtx_z, t, edep` |
| `truth_pf` | `<prefix>/pf` | double (Np, 25) | the 19 SBND columns, then `orig_trackid, gen_label, pf_node, cls, merged_into, nvis_q` |

Units cm / ns / GeV (the SBND attacher's). The set metadata carries `truth_file`, `truth_source` and the file's `meta` JSON. Schema: `WireCellTruth/Schema.h`; loader: `WireCellTruth/TruthTable.h`.

## The truth file

`truth_<ident>.npz`, one event per file, a stored zip without zip64 (the toolkit's cnpy reads no zip64): `truth_nu`, `truth_pf`, `sed_ids (N,2) int32 = trackid, orig_trackid` per SimEnergyDeposit in vector order (the depo id of Stage A index mode), `sed_trackq (T,2)` visible electrons per trackid, `anc_map (K,2)` dropped-descendant to kept-ancestor track ids, and the string tables `columns_nu`, `columns_pf`, `process_names`, `gen_labels`, `meta` as (K,L) uint8 ASCII. The `process` column is the index into `Truth::process_names()`; a file written with another list is remapped on load.

Design: `wcp-porting-img/pr_ml/docs/19_native_truth_plumbing.md` sec 3.4 (T4).
