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
  seen      the deposit lies inside the active volume (its face's wire extent) and the readout window, its main pixel (largest
            share) holds >= THR total true charge, and its channel is not bad at that tick
  free      (2.2) CS >= 0.5 in >= 1 seen view
  pieces    (2.5) per view, connected components (8-neighbour) of an object's own present pixels; a deposit joins the component
            of its main present entry (50f); 3-D pieces = 2-of-3 keys of the per-view components.  The extent of a piece is
            measured over its deposits SEEN in >= 2 views (G2), so a sub-hit deposit joined through a spread tail cannot inflate it
"""
import hashlib
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
    share, channel), each deposit's main entry, and `seen` with the bad channels applied.

    Speed-up pass (pr_ml doc 54 sec 6, byte-identical to the first build): narrow entry dtypes, the pixel key as a dense
    index over the view's touched channel range x padded slices, presence tables as dense arrays (direct gathers instead of
    searchsorted on the ring), the (pixel, instance) group table shared between clean-share calls with the same instances,
    and the ring test in chunks of entries.  Every sum is accumulated in the same entry order as before."""

    CHUNK = 1 << 22          # entries per ring-test chunk (bounds the transient memory; the sums are taken once over all entries)

    def __init__(self, ev, bad=None):
        D, geo = ev.D, ev.geo
        S = ev.projection_settings()
        self.n = ev.n
        self.tk = (D[:, 0] + S["toff_ns"]) / S["tick"] - S["tbin"]
        self.nslices = int(S["nticks"]) // 4
        self.NS = self.nslices + 2                                  # padded slice axis: index = slice + 1, ring slices -1 / nslices land in the pads
        self.inside = ev.inside
        self.sig_t = np.maximum(np.sqrt((D[:, 5] / S["speed"] / (S["tick"] / 1000.0)) ** 2 + (S["ext_long_us"] / (S["tick"] / 1000.0)) ** 2), 0.05)
        self.fp_list = [(an, fi, pi) for an, faces in geo.items() for fi, pl in faces.items() for pi in pl]
        self.fp_index = {k: i for i, k in enumerate(self.fp_list)}
        self.chan_of = {self.fp_index[k]: np.asarray(geo[k[0]][k[1]][k[2]]["chan"], np.int64) for k in self.fp_list}
        self.nwmax_of = np.array([len(self.chan_of[i]) for i in range(len(self.fp_list))], np.int64)
        self.NF = len(self.fp_list)
        # one flat channel table over all face planes: channel of (fp, wire) = chan_all[chan_off[fp] + wire]
        self.chan_off = np.r_[0, np.cumsum(self.nwmax_of)[:-1]].astype(np.int64)
        self.chan_all = np.concatenate([self.chan_of[i] for i in range(self.NF)]) if self.NF else np.zeros(0, np.int64)
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
                ok = (W > 0) & (WI >= 0) & (WI < len(g["chan"])) & (SI >= 0) & (SI < self.nslices)     # inside the readout window
                r_, k_ = np.nonzero(ok)
                if not len(r_):         # this face's deposits all outside the plane's wires / window: no entry, so a plane no deposit
                    del WI, SI, W, ok, r_, k_   # reaches stays None (pr_ml doc 58 sec 11.1; a zero-row entry failed _totals' ch.min())
                    continue
                wi = WI[r_, k_]
                fp = self.fp_index[(an, fi, pi)]
                ent.append((sel[r_].astype(np.int32), np.full(len(r_), fp, np.int16), wi.astype(np.int32), SI[r_, k_].astype(np.int32),
                            W[r_, k_], self.chan_of[fp][wi].astype(np.int32)))
                del WI, SI, W, ok, r_, k_, wi
            if not ent:
                self.ENT[pi] = None
                continue
            di, fpi, w, sl, wt, ch = (np.concatenate([e[x] for e in ent]) for x in range(6))
            del ent
            self.ENT[pi] = dict(di=di, fpi=fpi, w=w, sl=sl, wt=wt, ch=ch)
        self.q = ev.q
        self._grp = {}                                              # the (pixel, instance) group tables of the latest instance array
        self.itype = np.int32                                       # instance ids in the dense tables (checked against 2^31 in _groups)
        self._totals()
        self.set_bad(bad if bad is not None else ev.bad)

    # ------------------------------------------------------------------------------------------ pixel index
    def pixidx(self, pi, ch, sl, ds=0):
        """dense pixel index of (channel, slice + ds) in view pi: (channel - chmin) * NS + slice + ds + 1 (int64)"""
        return (ch.astype(np.int64) - self.chmin[pi]) * self.NS + (sl.astype(np.int64) + ds + 1)

    def _totals(self):
        """per view: the total true charge of each deposit's main pixel (largest share; first entry among equal shares), and the
        main entry itself.  A deposit's entries are contiguous in the entry arrays (one spread row), which the reduceat uses."""
        self.main = {}            # per view: (deposit indices with an entry, ascending; their main entry index)
        self.main_tot = {}        # per view: the main pixel's total true charge, aligned with main[pi][0]
        self.chmin = {}; self.nch = {}
        self.has_entry = np.zeros((self.n, 3), bool)
        for pi in range(3):
            E = self.ENT[pi]
            if E is None:
                self.main[pi] = (np.zeros(0, np.int64), np.zeros(0, np.int64)); self.main_tot[pi] = np.zeros(0)
                self.chmin[pi] = 0; self.nch[pi] = 0
                continue
            di, wt, ch, sl = E["di"], E["wt"], E["ch"], E["sl"]
            self.has_entry[di, pi] = True
            self.chmin[pi] = int(ch.min()); self.nch[pi] = int(ch.max()) - self.chmin[pi] + 1
            pix = self.pixidx(pi, ch, sl)
            ptot = np.bincount(pix, weights=wt * self.q[di], minlength=self.nch[pi] * self.NS)
            starts = np.flatnonzero(np.r_[True, di[1:] != di[:-1]])
            mx = np.maximum.reduceat(wt, starts)
            gid = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, len(di)]))
            idx = np.where(wt == mx[gid], np.arange(len(di)), len(di))
            me = np.minimum.reduceat(idx, starts)
            md = di[me].astype(np.int64)
            o = np.argsort(md, kind="stable")
            md, me = md[o], me[o]
            self.main[pi] = (md, me); self.main_tot[pi] = ptot[pix[me]]
            del pix, ptot, gid, idx

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
            ok = (self.main_tot[pi] >= THR) & self.inside[md]          # a deposit outside the active volume is never seen
            self.seen_nobad[md, pi] = ok
            isbad = np.zeros(len(me), bool)
            if self.bad:
                ch_m = E["ch"][me].astype(np.int64); tick_m = 4 * E["sl"][me].astype(np.int64) + 2
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
    def _groups(self, inst, pi):
        """the (pixel, instance) groups of view pi for the deposits with inst >= 0: pixel index, instance, own charge sum
        (entries accumulated in entry order, as the first build); the three views are cached for the latest instance array"""
        h = hashlib.blake2b(np.ascontiguousarray(inst).tobytes(), digest_size=16).digest()
        if self._grp.get("hash") != h:
            self._grp = {"hash": h}
        if pi in self._grp:
            return self._grp[pi]
        E = self.ENT[pi]
        di = E["di"]; ok_ = inst[di]; has = ok_ >= 0
        if not has.any():
            G = (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0))
        else:
            uinst = np.unique(inst[inst >= 0])                          # compact ids in instance order (deposit-level unique, cheap)
            if int(uinst[-1]) >= np.iinfo(np.int32).max:
                self.itype = np.int64
            look = np.full(int(uinst[-1]) + 1, -1, np.int64); look[uinst] = np.arange(len(uinst))
            ic = look[ok_[has]]
            pix = self.pixidx(pi, E["ch"][has], E["sl"][has])
            comb = pix * len(uinst) + ic
            del pix, ic
            uk, inv = np.unique(comb, return_inverse=True)
            del comb
            gq = np.bincount(inv.ravel(), weights=E["wt"][has] * self.q[di[has]])
            G = (uk // len(uinst), uinst[uk % len(uinst)], gq)
            del uk, inv
        self._grp[pi] = G
        return G

    def _presence(self, inst, pi, imask=None):
        """dense tables over the view's channel pixels: P1 = the present instance with the largest own charge (-1 none), P2 = the
        second (-1 none); ties by instance id ascending (the first build's lexsort order)"""
        gpix, ginst, gq = self._groups(inst, pi)
        pres = gq >= THR
        if imask is not None:
            pres &= imask[ginst]
        pc, po, pq_ = gpix[pres], ginst[pres], gq[pres]
        ncell = self.nch[pi] * self.NS
        P1 = np.full(ncell, -1, self.itype); P2 = np.full(ncell, -1, self.itype)
        if len(pc):
            o3 = np.lexsort((-pq_, pc)); pc, po = pc[o3], po[o3]
            upc, st_, cnt = np.unique(pc, return_index=True, return_counts=True)
            P1[upc] = po[st_]
            two = cnt >= 2
            P2[upc[two]] = po[st_[two] + 1]
        return P1, P2

    def _ring_table(self, pi, P1, P2):
        """the ring summary on the wire grid of the view's face planes with entries (one pad wire on each side of every face,
        the padded slice axis): R = -1 when no present instance lies in the 3x3 ring of the cell, the instance when exactly one
        does (over the P1 / P2 of the nine cells), -2 when two or more do.  A deposit's entry is clean iff R in (-1, its instance),
        which is the first build's test (unclean iff some ring cell has P1 or P2 outside {-1, k}).  Built by scattering the present
        channel pixels (a present channel is present at every wire reading it) onto their 3x3 neighbours.  -> (R, woff) with
        woff[fp] the padded wire offset (-1 for a face plane without entries)"""
        E = self.ENT[pi]
        fps = np.unique(E["fpi"]).astype(np.int64)
        woff = np.full(self.NF, -1, np.int64)
        nrow = 0
        for fp in fps.tolist():
            woff[fp] = nrow; nrow += int(self.nwmax_of[fp]) + 2
        NS = self.NS; chmin = self.chmin[pi]
        pc = np.flatnonzero(P1 >= 0)                                 # present channel cells
        chv = pc // NS + chmin; sv = pc % NS
        tg, vl = [], []
        for fp in fps.tolist():
            chans = self.chan_of[fp]
            order = np.argsort(chans, kind="stable"); cs = chans[order]
            lo = np.searchsorted(cs, chv); hi = np.searchsorted(cs, chv, side="right")
            cnt = hi - lo
            m = np.flatnonzero(cnt > 0)
            if not len(m):
                continue
            rep = np.repeat(m, cnt[m])                               # present cell index per (cell, wire) pair
            within = np.arange(len(rep)) - np.repeat(np.cumsum(cnt[m]) - cnt[m], cnt[m])
            wires = order[lo[rep] + within]
            wc = (woff[fp] + 1 + wires) * NS + sv[rep]               # the present wire cells of this face plane
            v1 = P1[pc[rep]]; v2 = P2[pc[rep]]
            for dw in (-1, 0, 1):
                for ds in (-1, 0, 1):
                    t = wc + dw * NS + ds
                    tg.append(t); vl.append(v1)
                    k2 = v2 >= 0
                    if k2.any():
                        tg.append(t[k2]); vl.append(v2[k2])
        R = np.full(nrow * NS, -1, self.itype)
        if tg:
            T_ = np.concatenate(tg); V_ = np.concatenate(vl)
            o = np.argsort(T_, kind="stable"); T_, V_ = T_[o], V_[o]
            st = np.flatnonzero(np.r_[True, T_[1:] != T_[:-1]])
            vmin = np.minimum.reduceat(V_, st); vmax = np.maximum.reduceat(V_, st)
            R[T_[st]] = np.where(vmin == vmax, vmin, -2).astype(self.itype)
            del T_, V_, o, st, vmin, vmax
        return R, woff

    def clean_share(self, inst, test, qmask=None):
        """2.2 for every deposit with test True: per view the clean share CS (n, 3) for the deposit's own instance inst[d]
        (inst < 0: no instance, never present, tested against every present instance), and cause (n, 3): a present instance
        other than the deposit's in its main pixel's ring (-1 none; the first found in the first build's ring order).  qmask (n,)
        multiplies the charge that can be present (e.g. shower instances only -> 'shower-ringed' shares); it must be constant
        within an instance."""
        n = self.n
        CS = np.zeros((n, 3)); cause = np.full((n, 3), -1, np.int64)
        imask = None
        if qmask is not None:
            qm = np.asarray(qmask, np.float64); hasi = inst >= 0
            nin = int(inst[hasi].max()) + 1 if hasi.any() else 0
            lo = np.full(nin, np.inf); hi = np.full(nin, -np.inf)
            np.minimum.at(lo, inst[hasi], qm[hasi]); np.maximum.at(hi, inst[hasi], qm[hasi])
            if not np.all((lo == hi) | ~np.isfinite(lo)):
                raise ValueError("clean_share: qmask must be constant within an instance")
            if not np.all(np.isin(hi[np.isfinite(hi)], (0.0, 1.0))):
                raise ValueError("clean_share: qmask must be 0 or 1")
            imask = hi > 0.5
        for pi in range(3):
            E = self.ENT[pi]
            if E is None:
                continue
            di, fpi, w, sl, wt = E["di"], E["fpi"], E["w"], E["sl"], E["wt"]
            T = np.flatnonzero(test[di])
            if not len(T):
                continue
            P1, P2 = self._presence(inst, pi, imask)
            if not (P1 >= 0).any():                              # nothing present in this view: every pixel is clean
                CS[:, pi] = np.bincount(di[T], weights=wt[T], minlength=n)
                del P1, P2
                continue
            # the clean test of every tested entry from the ring table
            R, woff = self._ring_table(pi, P1, P2)
            clean = np.zeros(len(T), bool)
            for c0 in range(0, len(T), self.CHUNK):
                Tc = T[c0:c0 + self.CHUNK]
                cell = (woff[fpi[Tc].astype(np.int64)] + 1 + w[Tc].astype(np.int64)) * self.NS + sl[Tc].astype(np.int64) + 1
                r = R[cell]
                clean[c0:c0 + len(Tc)] = (r == -1) | (r == inst[di[Tc]])
                del cell, r
            del R
            CS[:, pi] = np.bincount(di[T], weights=wt[T] * clean, minlength=n)
            del clean
            # cause at the main entry of each tested deposit: the first other present instance in the ring, in the first build's order
            md, me = self.main[pi]
            mm = test[md]
            dE, fE, wE, sE = md[mm], fpi[me[mm]].astype(np.int64), w[me[mm]].astype(np.int64), sl[me[mm]].astype(np.int64)
            kE = inst[dE]
            nwm = self.nwmax_of[fE]; off = self.chan_off[fE]
            chmin, chmax = self.chmin[pi], self.chmin[pi] + self.nch[pi] - 1
            oth = np.full(len(dE), -1, np.int64)
            for dw in (-1, 0, 1):
                wn = wE + dw
                valid = (wn >= 0) & (wn < nwm)
                wn = np.clip(wn, 0, nwm - 1)
                chn = self.chan_all[off + wn]
                valid &= (chn >= chmin) & (chn <= chmax)
                base = (np.clip(chn, chmin, chmax) - chmin) * self.NS + sE + 1
                for ds in (-1, 0, 1):
                    cell = base + ds
                    p1 = P1[cell]
                    hit = (p1 >= 0) & valid
                    a1 = hit & (p1 != kE)
                    p2 = P2[cell]
                    a2 = hit & ~a1 & (p2 >= 0) & (p2 != kE)
                    upd = (oth < 0) & a1; oth[upd] = p1[upd]
                    upd = (oth < 0) & a2; oth[upd] = p2[upd]
            cause[dE, pi] = oth
            del P1, P2, oth
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
