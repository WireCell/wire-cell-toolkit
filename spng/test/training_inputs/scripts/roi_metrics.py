"""ROI-wise efficiency and purity for time-contiguous ROIs.

An ROI is a maximal run of "activated" pixels that is contiguous along the time
axis (the last dimension) within a single channel.  Truth ROIs come from
``labels == 1``, reco ROIs from ``y > threshold``.

    efficiency = # true ROIs containing >= 1 reco pixel  /  # true ROIs
    purity     = # reco ROIs containing >= 1 true pixel  /  # reco ROIs

Tensors are shaped ``(batch, channel, time)``; the batch dimension is handled
without any python-level loop, and extra middle dimensions are allowed
(anything between the first and last dim is treated as more channels).  A 2-D
``(channel, time)`` tensor is accepted as a single-sample batch.
"""

import torch

__all__ = [
    "roi_table",
    "roi_efficiency",
    "roi_purity",
    "roi_metrics",
    "threshold_scan",
]


def _as_batched(x):
    """Reshape ``(B, ..., T)`` (or ``(C, T)``) to ``(B, rows, T)``."""
    if x.dim() < 2:
        raise ValueError(f"need at least 2 dims (channel, time), got {tuple(x.shape)}")
    if x.dim() == 2:
        return x.unsqueeze(0)
    return x.reshape(x.shape[0], -1, x.shape[-1])


def roi_table(mask, other=None):
    """Enumerate the ROIs in ``mask`` and (optionally) tag the matched ones.

    Parameters
    ----------
    mask : bool tensor, ``(B, ..., T)``
        Activated pixels.  ROIs are runs contiguous along the last dim.
    other : bool tensor, same shape, optional
        The comparison mask.  An ROI is "matched" if it contains at least one
        pixel that is set in ``other``.

    Returns
    -------
    dict of 1-D tensors, one entry per ROI (ordered by sample, row, time):
        ``sample``  batch index
        ``row``     flattened channel index within the sample
        ``start``   first time bin (inclusive)
        ``length``  width in time bins
        ``matched`` bool, present only when ``other`` is given
    """
    m = _as_batched(mask)
    batch, rows, nt = m.shape
    flat = m.reshape(-1, nt)

    # A run starts where the mask turns on.  Padding each row with False on the
    # left keeps runs from being merged across the row boundary.
    prev = torch.cat([torch.zeros_like(flat[:, :1]), flat[:, :-1]], dim=1)
    starts = flat & ~prev

    # Row-major cumsum over start flags labels every activated pixel with the
    # 1-based id of the ROI it belongs to (0 for inactive pixels).
    ids = torch.cumsum(starts.reshape(-1).to(torch.int64), 0).reshape(flat.shape)
    ids = torch.where(flat, ids, torch.zeros_like(ids))
    nroi = int(starts.sum())

    dev = m.device
    row_idx = torch.arange(batch * rows, device=dev).unsqueeze(1).expand_as(flat)
    col_idx = torch.arange(nt, device=dev).unsqueeze(0).expand_as(flat)

    out = {
        "sample": torch.div(row_idx[starts], rows, rounding_mode="floor"),
        "row": row_idx[starts] % rows,
        "start": col_idx[starts],
        # bincount over pixel ids; drop the id-0 (inactive) bucket
        "length": torch.bincount(ids[flat], minlength=nroi + 1)[1:],
    }

    if other is not None:
        o = _as_batched(other).reshape(-1, nt)
        if o.shape != flat.shape:
            raise ValueError("mask and other must have the same shape")
        hit = torch.zeros(nroi + 1, dtype=torch.bool, device=dev)
        hit[ids[flat & o]] = True
        out["matched"] = hit[1:]

    return out


def _rate(matched, sample, nsample, per_sample):
    """matched-ROI fraction, either pooled or per batch entry (nan if empty)."""
    if per_sample:
        tot = torch.bincount(sample, minlength=nsample).to(torch.float64)
        hit = torch.bincount(sample[matched], minlength=nsample).to(torch.float64)
        return torch.where(tot > 0, hit / tot, torch.full_like(tot, float("nan")))
    tot = matched.numel()
    if tot == 0:
        return torch.tensor(float("nan"), dtype=torch.float64, device=matched.device)
    return matched.sum().to(torch.float64) / tot


def roi_efficiency(y, labels, threshold=0.5, per_sample=False):
    """Fraction of true ROIs containing at least one pixel with ``y > threshold``.

    Returns a scalar tensor, or a ``(batch,)`` tensor if ``per_sample=True``
    (nan for samples with no true ROIs).
    """
    t = roi_table(labels == 1, y > threshold)
    return _rate(t["matched"], t["sample"], _as_batched(y).shape[0], per_sample)


def roi_purity(y, labels, threshold=0.5, per_sample=False):
    """Fraction of reco ROIs (``y > threshold``) containing at least one true pixel.

    Returns a scalar tensor, or a ``(batch,)`` tensor if ``per_sample=True``
    (nan for samples with no reco ROIs).
    """
    t = roi_table(y > threshold, labels == 1)
    return _rate(t["matched"], t["sample"], _as_batched(y).shape[0], per_sample)


def roi_metrics(y, labels, threshold=0.5, per_sample=False):
    """Efficiency and purity together, with the underlying ROI counts.

    Returns a dict with ``efficiency``, ``purity``, ``n_true``, ``n_true_matched``,
    ``n_reco``, ``n_reco_matched`` (counts are per-sample when requested).
    """
    nsample = _as_batched(y).shape[0]
    true_t = roi_table(labels == 1, y > threshold)
    reco_t = roi_table(y > threshold, labels == 1)

    def counts(tab):
        if per_sample:
            n = torch.bincount(tab["sample"], minlength=nsample)
            k = torch.bincount(tab["sample"][tab["matched"]], minlength=nsample)
            return n, k
        return (
            torch.tensor(tab["matched"].numel()),
            tab["matched"].sum().to(torch.int64),
        )

    n_true, k_true = counts(true_t)
    n_reco, k_reco = counts(reco_t)
    return {
        "threshold": threshold,
        "efficiency": _rate(true_t["matched"], true_t["sample"], nsample, per_sample),
        "purity": _rate(reco_t["matched"], reco_t["sample"], nsample, per_sample),
        "n_true": n_true,
        "n_true_matched": k_true,
        "n_reco": n_reco,
        "n_reco_matched": k_reco,
        "true_t": true_t,
        "reco_t": reco_t,
    }


def threshold_scan(y, labels, thresholds):
    """Run :func:`roi_metrics` over a sequence of thresholds.

    The truth ROIs are threshold independent, so they are enumerated once.
    Returns a dict of 1-D tensors aligned with ``thresholds``.
    """
    thresholds = [float(t) for t in thresholds]
    true_mask = labels == 1
    true_tab = roi_table(true_mask)
    true_sample, n_true = true_tab["sample"], len(true_tab["start"])
    nsample = _as_batched(y).shape[0]

    eff, pur, n_reco = [], [], []
    for t in thresholds:
        reco_mask = y > t
        hit = roi_table(true_mask, reco_mask)["matched"]
        reco = roi_table(reco_mask, true_mask)
        eff.append(float(_rate(hit, true_sample, nsample, False)))
        pur.append(float(_rate(reco["matched"], reco["sample"], nsample, False)))
        n_reco.append(len(reco["start"]))
    return {
        "threshold": torch.tensor(thresholds),
        "efficiency": torch.tensor(eff),
        "purity": torch.tensor(pur),
        "n_reco": torch.tensor(n_reco),
        "n_true": n_true,
    }
