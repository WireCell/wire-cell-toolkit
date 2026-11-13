"""Pure-torch replacements for the two torch_scatter functions DV2Net uses (torch_scatter is not installed in the
in-chain environment).  Semantics of torch_scatter 2.1 for dim=0 and a 1-D index:
  scatter_mean  sum per group / max(count, 1)   (empty groups -> 0)
  scatter_max   max per group, empty groups -> 0; returns (values, None) (DV2Net only reads [0])
"""
import torch


def _expand(index, src):
    return index.view(-1, *([1] * (src.dim() - 1))).expand_as(src)


def scatter_sum(src, index, dim=0, dim_size=None):
    if dim != 0:
        raise ValueError("only dim=0 is supported")
    n = int(index.max()) + 1 if dim_size is None else int(dim_size)
    out = src.new_zeros((n,) + tuple(src.shape[1:]))
    return out.index_add_(0, index.long(), src)


def scatter_mean(src, index, dim=0, dim_size=None):
    out = scatter_sum(src, index, dim, dim_size)
    cnt = torch.zeros(len(out), dtype=src.dtype, device=src.device).index_add_(
        0, index.long(), torch.ones(len(index), dtype=src.dtype, device=src.device))
    cnt = cnt.clamp(min=1)
    return out / cnt.view(-1, *([1] * (out.dim() - 1)))


def scatter_max(src, index, dim=0, dim_size=None):
    if dim != 0:
        raise ValueError("only dim=0 is supported")
    n = int(index.max()) + 1 if dim_size is None else int(dim_size)
    out = src.new_zeros((n,) + tuple(src.shape[1:]))
    if len(src):
        out.scatter_reduce_(0, _expand(index.long(), src), src, reduce="amax", include_self=False)
    return out, None
