"""labels_v2.pixels -- the 2-D side of the rules (50b secs 2.2, 2.5): one projection per event, then presence, ring,
clean share, pieces.

Projection (50d / 50f / 50g, kept exactly): a drifted deposit's charge is spread over the pixels (readout channel x 4-tick
slice) of each view of its face with erf shares: wire offsets -3..3 around rint(s), slice offsets -2..2 around floor(tk/4),
widths = the deposit's own drift sigmas (DepoFileSink row, diffusion applied) plus the campaign SP widths (truth.json
real_assoc transverse per plane, extra_smear_long in time) in quadrature; shares below 1e-4 dropped, the rest renormalised
to 1 per deposit per view; entries off the wire range dropped after that.

  present   an instance's own true charge in a pixel >= THR (2,400 e, the hit threshold)
  ring      the 3x3 neighbourhood: wire +-1 on the deposit's own face / plane, slice +-1
  clean     for instance k: no pixel of the ring holds a present instance other than k
  CS        clean share per view = sum of the deposit's shares over clean pixels
  seen      the deposit's main pixel (largest share) holds >= THR total true charge AND its channel is not bad at that tick
  free      (2.2) CS >= 0.5 in >= 1 seen view
  pieces    (2.5) per view, connected components (8-neighbour) of an object's own present pixels; a deposit joins the component
            of its main present entry; 3-D pieces = 2-of-3 keys of the per-view components (50f)
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.special import erf

THR = 2400.0
SQ2 = np.sqrt(2.0)
SOFF = 4096                 # slice offset so that codes stay positive
SL = 1 << 22                # pixel code = channel * SL + slice + SOFF
WB, SB = 1 << 13, 1 << 14   # wire / slice ranges in the per-object pixel key of the pieces


def spread(s, tk, sw, st, nw_off=3, ns_off=2):
    """per deposit: wire index (n, K), slice (n, K), share (n, K) normalised to 1 (d50d_real.spread, unchanged)"""
    w0 = np.rint(s).astype(np.int64); c0 = np.floor(tk / 4.0).astype(np.int64)
    dw = np.arange(-nw_off, nw_off + 1); ds = np.arange(-ns_off, ns_off + 1)
    ww = w0[:, None] + dw[None, :]
    fw = 0.5 * (erf((ww + 0.5 - s[:, None]) / SQ2 / sw[:, None]) - erf((ww - 0.5 - s[:, None]) / SQ2 / sw[:, None]))
    ss = c0[:, None] + ds[None, :]
    ft = 0.5 * (erf((4 * ss + 4 - tk[:, None]) / SQ2 / st[:, None]) - erf((4 * ss - tk[:, None]) / SQ2 / st[:, None]))
    W = (ft[:, :, None] * fw[:, None, :]).reshape(len(s), -1)
    WI = np.broadcast_to(ww[:, None, :], (len(s), len(ds), len(dw))).reshape(len(s), -1)
    SI = np.broadcast_to(ss[:, :, None], (len(s), len(ds), len(dw))).reshape(len(s), -1)
    W = np.where(W > 1e-4, W, 0.0)
    W = W / np.maximum(W.sum(1, keepdims=True), 1e-30)
    return WI, SI, W


class Projection:
    """the entries of every placed deposit in the three views: per view arrays (deposit, face-plane index, wire, slice,
    share, channel), the per-pixel total true charge, each deposit's main entry, and `seen` with the bad channels applied"""

    def __init__(self, ev, bad=None):
        D, geo = ev.D, ev.geo
        S = ev.projection_settings()
        self.n = ev.n
        self.tk = (D[:, 0] + S["toff_ns"]) / S["tick"] - S["tbin"]
        self.sig_t = np.maximum(np.sqrt((D[:, 5] / S["speed"] / (S["tick"] / 1000.0)) ** 2 + (S["ext_long_us"] / (S["tick"] / 1000.0)) ** 2), 0.05)
        self.fp_list = [(an, fi, pi) for an, faces in geo.items() for fi, pl in faces.items() for pi in pl]
        self.fp_index = {k: i for i, k in enumerate(self.fp_list)}
        self.chan_of = {self.fp_index[k]: np.asarray(geo[k[0]][k[1]][k[2]]["chan"], np.int64) for k in self.fp_list}
        self.nwmax_of = np.array([len(self.chan_of[i]) for i in range(len(self.fp_list))], np.int64)
        self.NF = len(self.fp_list)
        self.ENT = {}
        for pi in range(3):
            ent = []
            for (an, fi, pj) in self.fp_list:
                if pj != pi:
                    continue
                g = geo[an][fi][pi]
                sel = np.flatnonzero(ev.placed & (ev.da == an) & (ev.df == fi))
                if not len(sel):
                    continue
                s = np.einsum("ij,j->i", D[sel, 2:5] - g["cen0"], g["pdir"]) / g["pitch"]
                sw = np.maximum(np.sqrt(D[sel, 6] ** 2 + S["ext_tr"][pi] ** 2) / g["pitch"], 0.02)
                WI, SI, W = spread(s, self.tk[sel], sw, self.sig_t[sel])
                ok = (W > 0) & (WI >= 0) & (WI < len(g["chan"]))
                r_, k_ = np.nonzero(ok)
                ent.append((sel[r_], np.full(len(r_), self.fp_index[(an, fi, pi)], np.int64), WI[r_, k_], SI[r_, k_], W[r_, k_]))
            if not ent:
                self.ENT[pi] = None
                continue
            di, fpi, w, sl, wt = (np.concatenate([e[x] for e in ent]) for x in range(5))
            ch = np.empty(len(w), np.int64)
            for u in np.unique(fpi):
                m = fpi == u; ch[m] = self.chan_of[int(u)][w[m]]
            self.ENT[pi] = dict(di=di, fpi=fpi, w=w, sl=sl, wt=wt, ch=ch, code=ch * SL + sl + SOFF)
        self.q = ev.q
        self._totals()
        self.set_bad(bad if bad is not None else ev.bad)

    def _totals(self):
        """per view: the pixel codes (sorted) and the total true charge per pixel; each deposit's main entry (largest share)"""
        self.pix = {}
        self.main = {}            # per view: (deposit indices with an entry, their main entry index)
        for pi in range(3):
            E = self.ENT[pi]
            if E is None:
                self.pix[pi] = (np.zeros(0, np.int64), np.zeros(0)); self.main[pi] = (np.zeros(0, np.int64), np.zeros(0, np.int64))
                continue
            pu, pinv = np.unique(E["code"], return_inverse=True)
            ptot = np.bincount(pinv.ravel(), weights=E["wt"] * self.q[E["di"]])
            self.pix[pi] = (pu, ptot)
            o4 = np.lexsort((-E["wt"], E["di"])); f4 = np.r_[True, E["di"][o4][1:] != E["di"][o4][:-1]]
            self.main[pi] = (E["di"][o4][f4], o4[f4])

    def set_bad(self, bad):
        """seen[d, v]: main pixel total >= THR and channel not bad at tick 4 * slice + 2 (the slice centre)"""
        self.bad = bad or {}
        self.seen = np.zeros((self.n, 3), bool)
        self.seen_nobad = np.zeros((self.n, 3), bool)
        nbad_hits = 0
        for pi in range(3):
            E = self.ENT[pi]
            if E is None:
                continue
            md, me = self.main[pi]
            pu, ptot = self.pix[pi]
            ix = np.searchsorted(pu, E["code"][me])
            tot = ptot[np.minimum(ix, len(pu) - 1)]
            ok = tot >= THR
            self.seen_nobad[md, pi] = ok
            isbad = np.zeros(len(me), bool)
            if self.bad:
                ch_m = E["ch"][me]; tick_m = 4 * E["sl"][me] + 2
                o = np.argsort(ch_m, kind="stable"); cs = ch_m[o]
                for c, rngs in self.bad.items():
                    a, b = np.searchsorted(cs, c), np.searchsorted(cs, c, side="right")
                    if b <= a:
                        continue
                    sel = o[a:b]; t = tick_m[sel]
                    hit = np.zeros(len(sel), bool)
                    for t0, t1 in rngs:
                        hit |= (t >= t0) & (t < t1)
                    isbad[sel] |= hit
                nbad_hits += int(isbad.sum())
            self.seen[md, pi] = ok & ~isbad
        self.n_bad_main_entries = nbad_hits

    # ------------------------------------------------------------------------------------------ presence and the ring
    def clean_share(self, inst, test, qmask=None):
        """2.2 for every deposit with test True: per view the clean share CS (n, 3) for the deposit's own instance inst[d]
        (inst < 0: no instance, never present, tested against every present instance), and cause (n, 3): a present instance
        other than the deposit's in its main pixel's ring (-1 none).  qmask (n,) multiplies the charge that can be present
        (e.g. shower instances only -> 'shower-ringed' shares)."""
        n = self.n
        q = self.q if qmask is None else self.q * qmask
        CS = np.zeros((n, 3)); cause = np.full((n, 3), -1, np.int64)
        for pi in range(3):
            E = self.ENT[pi]
            if E is None:
                continue
            di, fpi, w, sl, wt, code = E["di"], E["fpi"], E["w"], E["sl"], E["wt"], E["code"]
            ok_ = inst[di]
            has = ok_ >= 0
            key = np.lexsort((ok_[has], code[has]))
            cc, oo, qq = code[has][key], ok_[has][key], (wt[has] * q[di[has]])[key]
            grp = np.r_[True, (cc[1:] != cc[:-1]) | (oo[1:] != oo[:-1])] if len(cc) else np.zeros(0, bool)
            gi = np.cumsum(grp) - 1
            gq = np.bincount(gi, weights=qq) if len(cc) else np.zeros(0)
            gc, go = cc[grp], oo[grp]
            pres = gq >= THR
            pc, po, pq_ = gc[pres], go[pres], gq[pres]
            o3 = np.lexsort((-pq_, pc)); pc, po = pc[o3], po[o3]
            upc, st_, cnt = np.unique(pc, return_index=True, return_counts=True)
            P1 = po[st_] if len(po) else np.zeros(0, np.int64)
            P2 = np.where(cnt >= 2, po[np.minimum(st_ + 1, max(len(po) - 1, 0))], -(1 << 60)) if len(po) else np.zeros(0, np.int64)
            T = np.flatnonzero(test[di])
            if not len(T):
                continue
            if not len(upc):                                   # nothing present in this view: every pixel is clean
                CS[:, pi] = np.bincount(di[T], weights=wt[T], minlength=n)
                continue
            dE, fE, wE, sE, wtE, kE = di[T], fpi[T], w[T], sl[T], wt[T], inst[di[T]]
            nwm = self.nwmax_of[fE]
            oth = np.full(len(T), -1, np.int64)
            for dw in (-1, 0, 1):
                wn = wE + dw
                valid = (wn >= 0) & (wn < nwm)
                wn = np.clip(wn, 0, nwm - 1)
                chn = np.empty(len(T), np.int64)
                for u in np.unique(fE):
                    m = fE == u; chn[m] = self.chan_of[int(u)][wn[m]]
                for ds in (-1, 0, 1):
                    cn = chn * SL + sE + ds + SOFF
                    ix = np.minimum(np.searchsorted(upc, cn), len(upc) - 1)
                    hit = (upc[ix] == cn) & valid
                    a1 = hit & (P1[ix] != kE)
                    a2 = hit & ~a1 & (cnt[ix] >= 2) & (P2[ix] != kE)
                    upd = (oth < 0) & a1; oth[upd] = P1[ix][upd]
                    upd = (oth < 0) & a2; oth[upd] = P2[ix][upd]
            clean = oth < 0
            CS[:, pi] = np.bincount(dE, weights=wtE * clean, minlength=n)
            # cause at the main entry of each tested deposit
            md, me = self.main[pi]
            pos = np.full(len(di), -1, np.int64); pos[T] = np.arange(len(T))
            mm = pos[me] >= 0
            cause[md[mm], pi] = oth[pos[me][mm]]
        return CS, cause

    def free(self, CS, test=None):
        f = ((CS >= 0.5) & self.seen).any(1)
        return f if test is None else (f & test)

    # ------------------------------------------------------------------------------------------ pieces (2.5)
    def components(self, dob, only=None):
        """per view: connected components of each object's own present pixels (8-neighbour); deposit -> the component of its
        main present entry (-1 none).  dob (n,): object index per deposit (-1 none).  only (n, bool): presence from all the
        object's deposits, but the graph and the assignment only from deposits with only True (d50g_real.components)"""
        n = self.n
        comp = np.full((n, 3), -1, np.int64)
        for pi in range(3):
            E = self.ENT[pi]
            if E is None:
                continue
            m = dob[E["di"]] >= 0
            di, fpi, w, sl, wt, ch = (E[k][m] for k in ("di", "fpi", "w", "sl", "wt", "ch"))
            if not len(di):
                continue
            ob = dob[di]
            ckey = (ob * (1 << 24) + ch) * SB + (sl + SOFF)
            uk, inv = np.unique(ckey, return_inverse=True)
            oq = np.bincount(inv.ravel(), weights=wt * self.q[di])
            pres = oq[inv.ravel()] >= THR
            if only is not None:
                pres = pres & only[di]
            pkey = ((ob * self.NF + fpi) * WB + w) * SB + (sl + SOFF)
            P = np.unique(pkey[pres])
            if not len(P):
                continue
            rows, cols = [np.arange(len(P))], [np.arange(len(P))]
            for dw in (-1, 0, 1):
                for ds in (-1, 0, 1):
                    if dw == 0 and ds == 0:
                        continue
                    nk = P + dw * SB + ds
                    ix = np.minimum(np.searchsorted(P, nk), len(P) - 1)
                    hit = P[ix] == nk
                    rows.append(np.flatnonzero(hit)); cols.append(ix[hit])
            rr, cc = np.concatenate(rows), np.concatenate(cols)
            A = coo_matrix((np.ones(len(rr), np.int8), (rr, cc)), shape=(len(P), len(P))).tocsr()
            _, lab = connected_components(A, directed=False)
            e_ok = np.flatnonzero(pres)
            o4 = e_ok[np.lexsort((-wt[e_ok], di[e_ok]))]
            first = np.r_[True, di[o4][1:] != di[o4][:-1]]
            md, me = di[o4][first], o4[first]
            comp[md, pi] = lab.astype(np.int64)[np.searchsorted(P, pkey[me])]
        return comp


def pieces_2of3(dob, comp, q):
    """3-D pieces from the per-view components (50f's rule): keys (object, cU, cV, cW) join when they share components in
    two views; a key seen in one view joins the highest-charge multi-view key sharing its component, else the first
    single-view key sharing it.  -> piece index per deposit (-1 none); pieces are numbered 0..P-1 over the event"""
    n = len(dob)
    em = np.flatnonzero((dob >= 0) & (comp >= 0).any(1))
    piece = np.full(n, -1, np.int64)
    if not len(em):
        return piece
    keys, kinv = np.unique(np.column_stack([dob[em], comp[em]]), axis=0, return_inverse=True)
    kinv = kinv.ravel(); kq = np.bincount(kinv, weights=q[em]); nseen = (keys[:, 1:] >= 0).sum(1)
    K = len(keys); er, ec = [], []
    for va, vb in ((0, 1), (0, 2), (1, 2)):
        m = np.flatnonzero((keys[:, 1 + va] >= 0) & (keys[:, 1 + vb] >= 0))
        if len(m) < 2:
            continue
        g = keys[m][:, [0, 1 + va, 1 + vb]]
        o_ = m[np.lexsort((g[:, 2], g[:, 1], g[:, 0]))]
        gg = keys[o_][:, [0, 1 + va, 1 + vb]]
        same = (gg[1:] == gg[:-1]).all(1)
        er.append(o_[:-1][same]); ec.append(o_[1:][same])
    single = np.flatnonzero(nseen == 1); multi = np.flatnonzero(nseen >= 2)
    for v in range(3):
        sv = single[keys[single, 1 + v] >= 0]
        if not len(sv):
            continue
        mv = multi[keys[multi, 1 + v] >= 0]
        best = {}
        for kx in mv[np.argsort(-kq[mv], kind="stable")].tolist():
            best.setdefault((int(keys[kx, 0]), int(keys[kx, 1 + v])), kx)
        first = {}
        for kx in sv.tolist():
            t = (int(keys[kx, 0]), int(keys[kx, 1 + v]))
            if t in best:
                er.append(np.array([kx])); ec.append(np.array([best[t]]))
            elif t in first:
                er.append(np.array([kx])); ec.append(np.array([first[t]]))
            else:
                first[t] = kx
    rr, cc = (np.concatenate(er), np.concatenate(ec)) if er else (np.zeros(0, np.int64), np.zeros(0, np.int64))
    _, klab = connected_components(coo_matrix((np.ones(len(rr), np.int8), (rr, cc)), shape=(K, K)).tocsr(), directed=False)
    piece[em] = klab.astype(np.int64)[kinv]
    return piece
