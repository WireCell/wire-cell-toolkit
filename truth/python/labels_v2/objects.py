"""labels_v2.objects -- the truth tree read for the rules (50b secs 2.5, 2.6, 2.7 G3 / G4 / G6): which G4 rows root an
object, of which kind, and which deposits belong to it.

Kinds of a root row:
  TRACK   a track kind (mu, pi, K, p, Sigma, nuclei): its own instance (2.2 decides Track / Stub); orphan flag (G4)
  DELTA   an electron born on a track kind by eIoni / muIoni / hIoni / ionIoni / muPairProd / hPairProd, or a
          muMinusCaptureAtRest electron below michel_capture_ke_min (2.6 candidate; Q12 B)
  MICHEL  the electron of a stopped muon: mu+ Decay with wait >= 1.5 ns, or muMinusCaptureAtRest with KE >= 5 MeV
  EM      a G3 EM object: an e / gamma whose mother is not e / gamma and that is not born on a track (primaries, pi0
          photons, photons of neutrons / hadron inelastics / mu brems, decays in flight, radiological electrons)
  OTHER   any other row with deposits (should be empty)
An object = its root plus the subtree reached through e / gamma only (em_fold).  wait = start(e) - start(mu) - |mu end -
mu start| / c.  A DELTA / MICHEL whose object holds < THR of charge gives no object: its electron deposits fold into the
parent track and each of its photon daughters roots an EM object of its own (2.6, "its radiation is a 2.5 object").
The same fold is applied by rules.py to every candidate that fails the resolved test.
"""
import numpy as np
import truth_table as TT

TRK = {13, 211, 321, 2212, 3222, 3112}
ONTRACK = {"eIoni", "muIoni", "hIoni", "ionIoni", "muPairProd", "hPairProd", "muMinusCaptureAtRest"}
KIND = {"TRACK": 0, "DELTA": 1, "MICHEL": 2, "EM": 3, "OTHER": 4}
KIND_NAME = {v: k for k, v in KIND.items()}
C_LIGHT = 299.792458     # mm / ns
WAIT_NS = 1.5
MICHEL_CAPTURE_KE_MIN = 5.0     # MeV
RADIO = ("Radioactivation", "RadioactiveDecay", "RadioactiveDecayBase")


class Truth:
    """the truth table as arrays, with the per-row kinds and the deposit -> row map"""

    def __init__(self, ev, thr):
        pf = ev.truth["truth_pf"]; c = TT.pf_col
        self.pf = pf; self.nrow = len(pf)
        self.tt = pf[:, c("trackid")].astype(np.int64)
        self.pdg = pf[:, c("pdg")].astype(np.int64); self.apdg = np.abs(self.pdg)
        self.proc = pf[:, c("process")].astype(np.int64)
        PN = TT.PROCESS_NAMES
        self.pname = np.array([PN[p] if 0 <= p < len(PN) else "unknown" for p in self.proc])
        self.mot = pf[:, c("mother_trackid")].astype(np.int64)
        self.row = {int(t): i for i, t in enumerate(self.tt.tolist())}
        self.mrow = np.array([self.row.get(int(m), -1) for m in self.mot.tolist()], np.int64)
        self.ke = pf[:, c("KE")].astype(np.float64) * 1e3            # MeV
        self.S0 = pf[:, [c("start_x"), c("start_y"), c("start_z")]] * 10.0   # mm (true frame)
        self.S1 = pf[:, [c("end_x"), c("end_y"), c("end_z")]] * 10.0
        self.t_start = pf[:, c("start_t")].astype(np.float64)        # ns
        self.nu_row = pf[:, c("nu_row")].astype(np.int64)
        self.gen = pf[:, c("gen_label")].astype(np.int64)
        self.gen_labels = ev.truth.get("gen_labels", [])
        mp = np.where(self.mrow >= 0, self.pdg[np.maximum(self.mrow, 0)], 0)
        self.mpdg = mp; amp = np.abs(mp)
        self.is_e = self.apdg == 11; self.is_g = self.pdg == 22; self.emp = self.is_e | self.is_g
        self.trk = np.isin(self.apdg, list(TRK)) | (self.apdg > 1_000_000_000)
        self.m_em = np.where(self.mrow >= 0, self.emp[np.maximum(self.mrow, 0)], False)
        trackmom = np.isin(amp, list(TRK)) | (amp > 1_000_000_000)
        # decays of muons
        self.wait = np.full(self.nrow, np.nan)
        dec = self.is_e & (self.pname == "Decay") & (amp == 13)
        for i in np.flatnonzero(dec):
            jm = self.mrow[i]
            self.wait[i] = self.t_start[i] - self.t_start[jm] - np.linalg.norm(self.S1[jm] - self.S0[jm]) / C_LIGHT
        self.dif = dec & ((self.wait < WAIT_NS) | (mp == 13))          # mu- Decay is always in flight
        cap = self.is_e & (self.pname == "muMinusCaptureAtRest") & (amp == 13)
        self.michel = (dec & ~self.dif) | (cap & (self.ke >= MICHEL_CAPTURE_KE_MIN))
        self.delta = self.is_e & trackmom & np.isin(self.pname, list(ONTRACK)) & ~self.michel
        self.ontrack = self.delta | self.michel
        self.g3 = self.emp & ~self.m_em & ~self.ontrack
        self.orphan = (self.apdg > 1_000_000_000) & (self.mrow < 0) & (self.pname != "primary")
        self.primary = self.pname == "primary"
        # deposit -> G4 row: the original id when in the table, else the stored (folded) id
        pairs, pinv = np.unique(np.column_stack([ev.otid, ev.tid]), axis=0, return_inverse=True)
        gp = np.array([self.row.get(int(o), self.row.get(int(t), -1)) for o, t in pairs.tolist()], np.int64)
        self.grow = gp[pinv.ravel()]                                  # (n,) row per deposit, -1 unknown
        self.q_row = np.bincount(np.maximum(self.grow, 0), weights=ev.q * (self.grow >= 0), minlength=self.nrow)
        self.e_row = np.bincount(np.maximum(self.grow, 0), weights=ev.energy * (self.grow >= 0), minlength=self.nrow)
        self.thr = thr
        # the kind of every row as a root
        self.kind = np.full(self.nrow, KIND["OTHER"], np.int64)
        self.kind[self.trk] = KIND["TRACK"]
        self.kind[self.g3] = KIND["EM"]
        self.kind[self.delta] = KIND["DELTA"]
        self.kind[self.michel] = KIND["MICHEL"]
        self.fold = np.zeros(self.nrow, bool)          # rows whose own deposits go to the mother's root (failed candidates)
        # children lists (CSR over rows sorted by mother row) for the folds
        o = np.argsort(self.mrow, kind="stable"); self._ch_rows = o
        self._ch_start = np.searchsorted(self.mrow[o], np.arange(self.nrow + 1))
        self.stopper = self.trk | self.g3 | self.ontrack | ~self.emp   # a chain stops at any non-EM row too

    # ------------------------------------------------------------------------------------------ roots
    def roots(self):
        """root row of every row: walk up through e / gamma mothers until a stopper that is not folded; a folded row's own
        deposits belong to its mother's root.  -> root_of (nrow,), -1 for a row with no root (orphan EM chain)"""
        root_of = np.full(self.nrow, -2, np.int64)
        for i0 in range(self.nrow):
            if root_of[i0] != -2:
                continue
            chain, i = [], i0
            while True:
                if root_of[i] != -2:
                    r = root_of[i]; break
                chain.append(i)
                if self.stopper[i] and not self.fold[i]:
                    r = i; break
                if self.mrow[i] < 0 or len(chain) > 10000:
                    r = -1 if self.emp[i] or self.fold[i] else i; break
                i = self.mrow[i]
            for k in chain:
                root_of[k] = r
        return root_of

    def fold_candidate(self, i):
        """a DELTA / MICHEL row i fails (below THR or unresolved): its electron deposits and those of its electron
        descendants (reached through electrons) go to the mother's root; every photon reached that way roots an EM object"""
        self.fold[i] = True
        self.kind[i] = KIND["OTHER"]
        stack = [i]
        while stack:
            j = stack.pop()
            for k in self._ch_rows[self._ch_start[j]:self._ch_start[j + 1]].tolist():
                if self.is_e[k] and not self.stopper[k]:
                    self.fold[k] = True; stack.append(k)
                elif self.is_g[k] and not self.stopper[k]:
                    self.stopper[k] = True; self.kind[k] = KIND["EM"]; self.g3[k] = True; self.radiation_of_fold = True
        self.stopper[i] = False

    def object_charge(self, root_of):
        """charge and deposited energy per root row"""
        ok = root_of >= 0
        q = np.bincount(np.where(ok, root_of, 0), weights=self.q_row * ok, minlength=self.nrow)
        e = np.bincount(np.where(ok, root_of, 0), weights=self.e_row * ok, minlength=self.nrow)
        return q, e

    def origin(self, i):
        """G6 origin of root row i: the generator label when stored, else the chain (d50i_real.origin)"""
        if self.dif[i]:
            return "DIF"
        if self.pname[i] in RADIO:
            return "radioactive"
        if self.mrow[i] < 0:
            g = self.gen[i]
            if 0 <= g < len(self.gen_labels):
                return "primary:" + self.gen_labels[g]
            return "primary"
        m = self.mrow[i]
        if abs(self.pdg[m]) == 111:
            return "pi0"
        x, k = m, 0
        while x >= 0 and k < 50:
            p = self.apdg[x]
            if p == 2112:
                return "neutron"
            if p in TRK or p > 1_000_000_000:
                return "muon" if p == 13 else "hadron"
            x = self.mrow[x]; k += 1
        return "other"

    def interaction(self, i):
        """interaction id: the nu_row of the primary ancestor (-1 none)"""
        x, k = i, 0
        while x >= 0 and k < 10000:
            if self.mrow[x] < 0:
                return int(self.nu_row[x])
            x = self.mrow[x]; k += 1
        return -1


def build(ev, thr):
    """-> Truth with the provisional roots: every TRACK / DELTA / MICHEL / EM row with >= THR of object charge is a root; a
    DELTA / MICHEL object below THR is folded (2.6).  Sets T.root_of, T.inst_dep (root row per deposit, -1 none)"""
    T = Truth(ev, thr)
    root_of = T.roots()
    q, _ = T.object_charge(root_of)
    weak = np.flatnonzero((T.kind == KIND["DELTA"]) | (T.kind == KIND["MICHEL"]))
    weak = weak[q[weak] < thr]
    for i in weak.tolist():
        T.fold_candidate(i)
    if len(weak):
        root_of = T.roots()
    T.root_of = root_of
    T.n_weak_folded = int(len(weak))
    T.inst_dep = np.where(T.grow >= 0, root_of[np.maximum(T.grow, 0)], -1)
    # a deposit whose row has no root (an orphan EM chain, or an id not in the table) gets its own fallback root by stored id
    lost = T.inst_dep < 0
    if lost.any():
        T.inst_dep = T.inst_dep.copy()
        T.inst_dep[lost] = T.nrow + ev.tid[lost]     # pseudo rows >= nrow
    return T
