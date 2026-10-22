# BlobDepoLabelTensor

Joins a cluster with the drifted depos that made it and emits the per-(blob, G4 track) charge table of `Truth::depo_blob_assoc` as tensors, optionally together with the filled cluster.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobDepoLabelTensor` |
| Concrete class | `WireCell::Truth::BlobDepoLabelTensor` |
| Node category | join |
| Primary interface | `IClusterDepoTensor` |
| Input type(s) | `ICluster`, `IDepoSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

The association knobs of `BlobDepoLabel` (`time_offset`, `nsigma`, `speed`, `pindex`, `id_from_prior`, `extra_smear_long`, `extra_smear_tran`) plus:

| Parameter | Description |
| --- | --- |
| `datapath` | datapath pattern of the output tensors, `%d` = cluster ident (default `clusters/%d`) |
| `with_cluster` | include the filled cluster as TensorDM `cluster` tensors (default true) |
| `min_frac` | drop (blob, track) rows whose charge is below this fraction of the blob's total (default 0) |

## Output tensors

| datatype | shape | columns (also in `md["columns"]`) |
| --- | --- | --- |
| `blobtruth` | double (M, 6) | `blob_ident, face, slice_ident, trackid, q_bk, e_bk`; one row per (blob, track) with charge |
| `blobtruth_sum` | double (B, 7) | `blob_ident, face, slice_ident, q_true, q_top, tid_top, ntracks`; one row per blob in graph vertex order (blob idents repeat across slices and faces: join on the triple); `tid_top = -1` for a blob with no true charge; `md.ndepos` / `md.q_depos` = the depos in the anode's sensitive volumes and their charge sum |
| `cluster`, `clnodeset`, `cledgeset` | TensorDM | the filled cluster, when `with_cluster` |

`e_bk` is the key depo's `energy()` times the charge share and is meaningful only when the depo file carried the truth sidecar (`DepoFileSink.truth_sidecar`). Schema: `WireCellTruth/Schema.h`.
