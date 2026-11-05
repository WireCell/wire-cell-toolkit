"""labels_v2.writer -- <event>.labels_v2.npz: a deterministic zip (members sorted, fixed date, deflate 6; as
pr_ml/s0/convert_event.write_npz) with the per-point arrays, the instance table, the key points and a json meta."""
import io, json, os, zipfile
import numpy as np

DATE = (1980, 1, 1, 0, 0, 0)


def write_npz(fn, arrays):
    tmp = fn + ".part"
    os.makedirs(os.path.dirname(os.path.abspath(fn)), exist_ok=True)
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for k in sorted(arrays):
            b = io.BytesIO(); np.lib.format.write_array(b, np.ascontiguousarray(arrays[k]), allow_pickle=False)
            zi = zipfile.ZipInfo(k + ".npy", date_time=DATE); zi.external_attr = 0o644 << 16
            z.writestr(zi, b.getvalue(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    os.replace(tmp, fn)


def narrow(a):
    """the smallest integer type holding the array (as slim_event.narrow)"""
    a = np.asarray(a)
    if a.dtype.kind not in "iu" or not a.size:
        return a
    lo, hi = int(a.min()), int(a.max())
    for t in (np.int8, np.int16, np.int32):
        if np.iinfo(t).min <= lo and hi <= np.iinfo(t).max:
            return a.astype(t)
    return a.astype(np.int64)


def assemble(ev, R, pt, kps, kpt, meta):
    A = {}
    for k in ("cls7", "cls_obj", "inst", "piece_id", "piece_type", "overlap", "stretch", "other_reason", "lab_tid", "lab_orig_tid", "lab_cls", "lab_pf_node", "cl", "dep"):
        A["pt_" + k] = narrow(pt[k])
    for k in ("d_own_seg", "match_dist", "ghost_dist", "lab_d_traj"):
        A["pt_" + k] = np.asarray(pt[k], np.float32)
    A["pt_kp_type"] = narrow(kpt["kp_type"]); A["pt_kp_dist"] = kpt["kp_dist"]; A["pt_kp_index"] = narrow(kpt["kp_index"])
    A["pt_x"] = ev.pts["X"][:, 0].astype(np.float32); A["pt_y"] = ev.pts["X"][:, 1].astype(np.float32); A["pt_z"] = ev.pts["X"][:, 2].astype(np.float32)
    T = R.table
    for k, v in T.items():
        if k.startswith("inst_"):
            A[k] = narrow(v) if np.asarray(v).dtype.kind in "iu" else np.asarray(v, np.float32)
    A["inst_origin_names"] = np.frombuffer(json.dumps(T["origins"]).encode(), np.uint8)
    for k, v in kps.items():
        A[k] = narrow(v) if np.asarray(v).dtype.kind in "iu" else np.asarray(v, np.float32)
    A["meta"] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    return A
