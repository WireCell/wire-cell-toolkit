"""Batch collation for inference: collate / to_device vendored verbatim from dv2max/data.py (the training-label
fields are filled with the neutral values infer._inference_sample sets)."""
import numpy as np
import torch


def collate(samples):
    B = len(samples)
    nv = [len(s["vox_grid"]) for s in samples]
    off = np.concatenate([[0], np.cumsum(nv)])
    vc = np.concatenate([np.concatenate([np.full((n, 1), i, np.int64), s["vox_grid"].astype(np.int64)], 1)
                         for i, (s, n) in enumerate(zip(samples, nv))])
    shape = np.max(np.stack([s["grid_shape"] for s in samples]), 0)
    K = max(1, max(len(s["visible"]) for s in samples))
    vis = np.zeros((B, K, 3), np.float32)
    vmask = np.zeros((B, K), bool)
    for i, s in enumerate(samples):
        k = len(s["visible"])
        vis[i, :k] = s["visible"]
        vmask[i, :k] = True
    H = max(1, max(len(s["hard_neg"]) for s in samples)) if "hard_neg" in samples[0] else 1
    hn = np.zeros((B, H, 3), np.float32)
    hmask = np.zeros((B, H), bool)
    for i, s in enumerate(samples):
        k = len(s.get("hard_neg", ()))
        if k:
            hn[i, :k] = s["hard_neg"]
            hmask[i, :k] = True
    b = {"vox_coords": torch.from_numpy(vc), "spatial_shape": [int(v) for v in shape], "batch_size": B,
         "vox_feat": torch.from_numpy(np.concatenate([s["vox_feat"] for s in samples])),
         "vox_coord": torch.from_numpy(np.concatenate([s["vox_coord"] for s in samples])),
         "vox_cluster": torch.from_numpy(np.concatenate([s["vox_cluster"] for s in samples])),
         "pt_feat": torch.from_numpy(np.concatenate([s["pt_feat"] for s in samples])),
         "pt2vox": torch.from_numpy(np.concatenate([s["pt2vox"] + off[i] for i, s in enumerate(samples)])),
         "pt_xyz": torch.from_numpy(np.concatenate([s["pt_xyz"] for s in samples])),
         "ev_feat": torch.from_numpy(np.stack([s["ev_feat"] for s in samples])),
         "has_target": torch.tensor([bool(s["has_target"]) for s in samples]),
         "target": torch.from_numpy(np.stack([s["target"] for s in samples])),
         "vis": torch.from_numpy(vis), "vis_mask": torch.from_numpy(vmask),
         "hard_neg": torch.from_numpy(hn), "hard_neg_mask": torch.from_numpy(hmask),
         "loc_ok": torch.tensor([bool(s["loc_ok"]) for s in samples]),
         "host_cluster": torch.tensor([int(s["host_cluster"]) for s in samples], dtype=torch.long),
         "crop_label": torch.tensor([bool(s["crop_label"]) for s in samples]),
         "keys": [s["key"] for s in samples]}
    return b


def to_device(b, device):
    return {k: (v.to(device, non_blocking=True) if torch.is_tensor(v) else v) for k, v in b.items()}
