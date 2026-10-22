#!/usr/bin/env python3
"""truth_table.py -- the truth file (truth_<ident>.npz) of WireCellTruth, python side.

The file holds, per event, the tables of truth/inc/WireCellTruth/Schema.h:

  truth_nu     (Nn, 12) float64   one row per generator interaction (SBND attacher columns)
  truth_pf     (Np, 25) float64   one row per G4 particle, ALL kept (no KE cut; PFMerger applies pf_ke_min)
  sed_ids      (Ns, 2)  int32     per SimEnergyDeposit in vector order: trackid, orig_trackid
                                  (the depo id of Stage A index mode is the row number)
  sed_trackq   (Nt, 2)  float64   per trackid: visible electrons summed over its deposits
  anc_map      (Na, 2)  int32     dropped descendant trackid -> kept ancestor trackid (ParticleAncestryMap)
  sed_energy   (Ns,)    float32   energy (MeV) per SimEnergyDeposit in vector order (the depo energy of index mode)
  columns_nu, columns_pf, process_names, gen_labels, meta
                        (K, L) uint8 ASCII matrices, zero padded (the toolkit's cnpy reads no unicode)

Written as a plain STORED zip without zip64 (WCT's cnpy cannot parse the zip64 local headers
numpy.savez forces), one event per file.  Units follow the SBND attacher: cm, ns, GeV.
Readers index columns by name (columns_* / md["columns"]), never by position.

This module is imported by dump_truth_table.py (Stage A, inside LArSoft) and by the validation
scripts; it needs numpy only.  `python3 -I truth_table.py --selftest [out.npz]` writes a synthetic
file and reads it back.
"""
import io, json, sys, zipfile
import numpy as np

SCHEMA_VERSION = 3

NU_COLUMNS = ["nu_idx", "pdg", "ccnc", "mode", "int_type", "flavor",
              "E", "vtx_x", "vtx_y", "vtx_z", "t", "edep"]

PF_COLUMNS = [  # SBND attacher columns first and in order, then the native additions
    "nu_row", "trackid", "parent_trackid", "mother_trackid", "pdg", "process",
    "E", "KE",
    "start_x", "start_y", "start_z", "start_t",
    "end_x", "end_y", "end_z", "end_t",
    "start_px", "start_py", "start_pz",
    "orig_trackid", "gen_label", "pf_node", "cls", "merged_into", "nvis_q"]

# Index = code; mirrors Truth::process_names() in truth/src/Schema.cxx.  Append only.
PROCESS_NAMES = [
    "unknown", "primary",
    "eIoni", "muIoni", "hIoni", "ionIoni",
    "eBrem", "muBrems", "hBrems",
    "annihil", "conv", "compt", "phot", "Rayl",
    "Decay", "muMinusCaptureAtRest", "muMinusAtomicCapture",
    "hadElastic", "CoulombScat", "msc",
    "neutronInelastic", "protonInelastic", "pi+Inelastic", "pi-Inelastic",
    "kaon+Inelastic", "kaon-Inelastic", "kaon0LInelastic", "kaon0SInelastic",
    "lambdaInelastic", "sigma+Inelastic", "sigma-Inelastic",
    "anti_protonInelastic", "anti_neutronInelastic",
    "dInelastic", "tInelastic", "He3Inelastic", "alphaInelastic", "ionInelastic",
    "nCapture", "hFritiofCaptureAtRest", "hBertiniCaptureAtRest",
    "muPairProd", "hPairProd", "photonNuclear", "electronNuclear", "muonNuclear",
    "RadioactiveDecay", "RadioactiveDecayBase", "Radioactivation",
    "Transportation", "StepLimiter"]
PROCESS_UNKNOWN = 0
CLS = {"track": 0, "shower": 1, "michel": 2, "delta": 3, "le": 4, "unknown": -1}
TID_GHOST, TID_SHARED = -1, -2


def process_code(name):
    try:
        return PROCESS_NAMES.index(name)
    except ValueError:
        return PROCESS_UNKNOWN


def nu_col(name):
    return NU_COLUMNS.index(name)


def pf_col(name):
    return PF_COLUMNS.index(name)


def strings_to_u8(strings):
    """list of str -> (K, L) uint8 ASCII, zero padded; (0, 1) when empty"""
    if not strings:
        return np.zeros((0, 1), dtype=np.uint8)
    bs = [s.encode("ascii", "replace") for s in strings]
    L = max(1, max(len(b) for b in bs))
    out = np.zeros((len(bs), L), dtype=np.uint8)
    for i, b in enumerate(bs):
        out[i, :len(b)] = np.frombuffer(b, dtype=np.uint8)
    return out


def u8_to_strings(arr):
    return [bytes(row).split(b"\0", 1)[0].decode("ascii") for row in np.asarray(arr, dtype=np.uint8)]


def savez_nozip64(fn, **arrs):
    """numpy arrays into a STORED zip without zip64 headers (readable by WCT's cnpy)"""
    with zipfile.ZipFile(fn, "w", zipfile.ZIP_STORED, allowZip64=False) as z:
        for k, a in arrs.items():
            b = io.BytesIO()
            np.lib.format.write_array(b, np.ascontiguousarray(a), version=(1, 0))
            z.writestr(k + ".npy", b.getvalue())


def write_truth_npz(fn, nu, pf, sed_ids=None, sed_trackq=None, anc_map=None, gen_labels=(), meta=None, sed_energy=None):
    nu = np.asarray(nu, dtype=np.float64).reshape(-1, len(NU_COLUMNS))
    pf = np.asarray(pf, dtype=np.float64).reshape(-1, len(PF_COLUMNS))
    sed_ids = np.zeros((0, 2), np.int32) if sed_ids is None else np.asarray(sed_ids, dtype=np.int32).reshape(-1, 2)
    sed_trackq = np.zeros((0, 2), np.float64) if sed_trackq is None else np.asarray(sed_trackq, dtype=np.float64).reshape(-1, 2)
    anc_map = np.zeros((0, 2), np.int32) if anc_map is None else np.asarray(anc_map, dtype=np.int32).reshape(-1, 2)
    sed_energy = np.zeros(0, np.float32) if sed_energy is None else np.asarray(sed_energy, dtype=np.float32).reshape(-1)
    meta = dict(meta or {})
    meta.setdefault("schema_version", SCHEMA_VERSION)
    savez_nozip64(fn, truth_nu=nu, truth_pf=pf, sed_ids=sed_ids, sed_trackq=sed_trackq, anc_map=anc_map, sed_energy=sed_energy,
                  columns_nu=strings_to_u8(NU_COLUMNS), columns_pf=strings_to_u8(PF_COLUMNS),
                  process_names=strings_to_u8(PROCESS_NAMES), gen_labels=strings_to_u8(list(gen_labels)),
                  meta=strings_to_u8([json.dumps(meta, sort_keys=True)]))


def read_truth_npz(fn):
    """-> dict with numpy tables and decoded string lists; meta as a dict"""
    f = np.load(fn)
    out = {k: f[k] for k in f.files}
    for k in ("columns_nu", "columns_pf", "process_names", "gen_labels", "meta"):
        if k in out:
            out[k] = u8_to_strings(out[k])
    out["meta"] = json.loads(out["meta"][0]) if out.get("meta") else {}
    assert out["columns_nu"] == NU_COLUMNS, out["columns_nu"]
    assert out["columns_pf"] == PF_COLUMNS, out["columns_pf"]
    return out


def selftest(fn):
    nu = np.zeros((1, len(NU_COLUMNS)))
    nu[0, nu_col("nu_idx")] = 0; nu[0, nu_col("pdg")] = 14; nu[0, nu_col("E")] = 1.5
    pf = np.zeros((2, len(PF_COLUMNS)))
    pf[0, pf_col("trackid")] = 1; pf[0, pf_col("pdg")] = 13; pf[0, pf_col("process")] = process_code("primary")
    pf[1, pf_col("trackid")] = 2; pf[1, pf_col("pdg")] = 11; pf[1, pf_col("process")] = process_code("Decay")
    pf[1, pf_col("mother_trackid")] = 1; pf[1, pf_col("parent_trackid")] = 1
    pf[:, pf_col("pf_node")] = -1; pf[:, pf_col("cls")] = CLS["unknown"]; pf[:, pf_col("merged_into")] = -1
    write_truth_npz(fn, nu, pf, sed_ids=[(1, 1), (2, 7), (2, 2)], sed_trackq=[(1, 1e5), (2, 3e3)],
                    gen_labels=["generator"], meta={"detector": "selftest", "id_mode": "index"})
    r = read_truth_npz(fn)
    assert r["truth_pf"].shape == (2, len(PF_COLUMNS)) and r["sed_ids"].tolist() == [[1, 1], [2, 7], [2, 2]]
    assert r["process_names"] == PROCESS_NAMES and r["gen_labels"] == ["generator"]
    assert r["meta"]["detector"] == "selftest" and r["meta"]["schema_version"] == SCHEMA_VERSION
    print("selftest ok:", fn, {k: getattr(v, "shape", len(v)) for k, v in r.items() if k != "meta"})


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--selftest":
        selftest(sys.argv[2] if len(sys.argv) > 2 else "truth_selftest.npz")
    else:
        print(__doc__)
