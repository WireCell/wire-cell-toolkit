'''
Analysis of the SNB TA study truth tables (TATruth output, .npz files from
trig/cfg/snb-ta-study.jsonnet).

Definitions (see trig/docs/technote):

- A TA is "signal" if signal (SNB) sources give at least `frac` of its
  matched charge.  A TA with any SNB charge is "SNB-related": it exists (or
  differs) only because of the boosted SNB rate, so its rate is divided by
  the boost.  "Background" TAs have no SNB charge.
- An SNB interaction is "visible" if its charge reaches the collection
  plane inside the readout window (t_collection more than `margin` from
  either end).  It is "found" if it gives at least `frac` of the charge of a
  selected TA.  Efficiency = found / visible.
- The boost is the simulated SNB rate over the physical reference rate per
  tonne times the simulated argon mass.
- Rates are per simulated unit (an HD APA or a VD CRU).
- Data rate: each selected TA seeds a block of `block` around its middle
  over all channels of the unit.  Overlapping blocks merge.  Background
  blocks give the steady state, SNB-related blocks (unboosted) add during a
  burst (sparse approximation).

Selection methods:

- "charge" :: a TA on the seed plane with adc_integral >= cut.
- "sum" :: a seed TA (adc_integral >= seed) whose neighborhood sum
  (adc_integral of the TAs on the same plane and face within `radius` in
  space, ie radius/pitch channels and radius/drift_speed in time, with
  adc_integral >= neighbor) is >= cut.

Requires only numpy.
'''

import json

import numpy

# WCT units
ns = 1.0
us = 1000.0 * ns
ms = 1000.0 * us
s = 1000.0 * ms
mm = 1.0
cm = 10.0 * mm


def load(paths):
    '''Return (tas, interactions, links, frames) concatenated over files.

    Tables are dicts of column name -> array, frames a list of the tensor set
    metadata.  Frame numbers are made unique over files.  The interactions
    get a "file" column, the index of their file in paths.'''
    tables = {"tas": [], "interactions": [], "links": []}
    frames = []
    offset = 0
    for ifile, path in enumerate(paths):
        f = numpy.load(path)
        idents = sorted({int(k.split("_")[1]) for k in f.keys() if k.startswith("tensorset_")})
        maxframe = -1
        for ident in idents:
            md = json.loads(f[f"tensorset_{ident}_metadata.json"])
            md["frame"] += offset
            maxframe = max(maxframe, md["frame"])
            frames.append(md)
            for ind in range(3):
                tmd = json.loads(f[f"tensor_{ident}_{ind}_metadata.json"])
                arr = f[f"tensor_{ident}_{ind}_array"]
                if arr.size == 0:
                    arr = numpy.zeros((0, len(tmd["columns"])))
                arr = arr.reshape(-1, len(tmd["columns"])).copy()
                arr[:, tmd["columns"].index("frame")] += offset
                cols = tmd["columns"]
                if tmd["name"] == "interactions":
                    cols = cols + ["file"]
                    arr = numpy.hstack([arr, numpy.full((len(arr), 1), ifile)])
                tables[tmd["name"]].append((cols, arr))
        offset = maxframe + 1

    def concat(parts):
        cols = parts[0][0]
        arr = numpy.concatenate([a for _, a in parts])
        return {c: arr[:, i] for i, c in enumerate(cols)}

    return concat(tables["tas"]), concat(tables["interactions"]), concat(tables["links"]), frames


def union_length(a, b):
    '''Total length of the union of the intervals [a, b).'''
    if len(a) == 0:
        return 0.0
    order = numpy.argsort(a, kind="stable")
    a, b = a[order], b[order]
    bmax = numpy.maximum.accumulate(b)
    new = numpy.ones(len(a), dtype=bool)
    new[1:] = a[1:] > bmax[:-1]
    starts = numpy.flatnonzero(new)
    ends = numpy.append(starts[1:], len(a)) - 1
    return float(numpy.sum(bmax[ends] - a[starts]))


class Study:
    """The tables of one run set (one detector configuration)."""

    def __init__(self, paths, frac=0.5, margin_us=50.0, stress=("nCapture",)):
        self.paths = list(paths)
        self.tas, self.ints, self.links, self.frames = load(self.paths)
        md0 = self.frames[0]
        self.meta = md0["metadata"]
        self.detector = self.meta.get("detector", "hd")
        self.sources = md0["sources"]
        self.source_names = [src["name"] for src in self.sources]
        self.nsources = len(self.source_names)
        self.nchannels = md0["nchannels"]
        self.mass_kg = md0.get("mass_kg", float("nan"))
        self.pitch = md0.get("pitch", [4.7 * mm] * 3)
        self.drift_speed = md0.get("drift_speed", 1.6 * mm / us)
        self.tick = self.meta["tick"]
        self.frac = frac
        self.stress = set(stress)

        # boosts: SNB from the reference rate, others from the metadata
        if "snb_physical_rate_hz_per_t" in self.meta:
            self.snb_physical_rate_hz = self.meta["snb_physical_rate_hz_per_t"] * self.mass_kg / 1000.0
        else:  # v1
            self.snb_physical_rate_hz = self.meta["snb_physical_rate_hz"]
        self.snb_rate_hz = self.meta["snb_rate_hz"]
        self.boost = self.snb_rate_hz / self.snb_physical_rate_hz
        given = self.meta.get("source_boosts", {})
        self.source_boost = numpy.ones(self.nsources)
        for i, src in enumerate(self.sources):
            if src.get("signal"):
                self.source_boost[i] = self.boost
            elif src["name"] in given:
                self.source_boost[i] = given[src["name"]]
        self.source_stress = numpy.array([n in self.stress for n in self.source_names], dtype=bool)

        self.nframes = len(self.frames)
        nfr = max(fr["frame"] for fr in self.frames) + 1
        self.lo = numpy.zeros(nfr)
        self.hi = numpy.zeros(nfr)
        for fr in self.frames:
            self.lo[fr["frame"]] = fr["time_begin"]
            self.hi[fr["frame"]] = fr["time_end"]
        self.readout = float(numpy.median(self.hi - self.lo))
        self.live = sum(fr["time_end"] - fr["time_begin"] for fr in self.frames) / s

        t = self.tas
        with numpy.errstate(invalid="ignore", divide="ignore"):
            self.sfrac = numpy.where(t["q_total"] > 0, t["q_signal"] / t["q_total"], 0.0)
        self.ta_signal = self.sfrac >= frac
        self.ta_snb = t["q_signal"] > 0
        self.ta_source = t["dominant_source"].astype(int)  # -1: noise
        self.ta_dom = t["dominant_id"].astype(numpy.int64)
        known = self.ta_source >= 0
        src = numpy.where(known, self.ta_source, 0)
        # the boost that made each TA: SNB if it has SNB charge, else that of
        # its dominant source
        self.ta_boost = numpy.where(self.ta_snb, self.boost,
                                    numpy.where(known, self.source_boost[src], 1.0))
        self.ta_stress = ~self.ta_snb & known & self.source_stress[src]
        self.ta_boosted = self.ta_boost > 1

        taframe = t["frame"].astype(numpy.int64)
        if numpy.any(numpy.diff(taframe) < 0):
            raise ValueError("TA rows are not in frame order")
        first = numpy.searchsorted(taframe, numpy.arange(nfr))
        L = self.links
        lframe = L["frame"].astype(numpy.int64)
        self.link_ta = first[lframe] + L["ta"].astype(numpy.int64)

        ints = self.ints
        ikey = ints["frame"].astype(numpy.int64) * 2**32 + ints["id"].astype(numpy.int64)
        iorder = numpy.argsort(ikey)
        lkey = lframe * 2**32 + L["id"].astype(numpy.int64)
        pos = numpy.clip(numpy.searchsorted(ikey[iorder], lkey), 0, max(len(ikey) - 1, 0))
        ok = ikey[iorder][pos] == lkey if len(ikey) else numpy.zeros(len(lkey), dtype=bool)
        self.link_int = numpy.where(ok, iorder[pos], -1)
        self.link_dominant = ok & (L["q"] >= frac * t["q_total"][self.link_ta])

        margin = margin_us * us
        ifr = ints["frame"].astype(numpy.int64)
        tc = ints["t_collection"]
        inwin = (tc >= self.lo[ifr] + margin) & (tc <= self.hi[ifr] - margin)
        self.int_source = ints["source"].astype(int)
        # signal interactions arriving in the readout, visible or not
        self.inwin = (ints["signal"] == 1) & inwin
        self.visible = self.inwin & (ints["q_collection"] > 0)
        # interactions are generated uniformly over [-pre, readout): the
        # number arriving in the readout is this fraction of all
        self.gen_span = self.readout - float(numpy.min(ints["time"])) if len(ints["time"]) else self.readout
        self.window_fraction = self.readout / self.gen_span

    # selections

    def plane_mask(self, plane):
        return numpy.ones(len(self.tas["plane"]), dtype=bool) if plane is None else self.tas["plane"] == plane

    def charge_selection(self, cut, plane=2, min_ninputs=1):
        t = self.tas
        return self.plane_mask(plane) & (t["adc_integral"] >= cut) & (t["ninputs"] >= min_ninputs)

    def neighborhood_sum(self, plane=2, radius=30 * cm, neighbor=0.0):
        """The adc_integral sum around each TA on the plane (NaN elsewhere).

        A neighbor TA made by a boosted source (SNB or a boosted background)
        counts only if it is from the same interaction as the TA: other
        boosted interactions are absent at physical rates.  Physical rate
        sources (39Ar, ...) and noise always count."""
        t = self.tas
        out = numpy.full(len(t["plane"]), numpy.nan)
        m = numpy.flatnonzero(self.plane_mask(plane))
        if len(m) == 0:
            return out
        dch = radius / self.pitch[plane if plane is not None else 2]
        dt = radius / self.drift_speed
        key = t["frame"][m].astype(numpy.int64) * 4 + t["face"][m].astype(numpy.int64)
        order = numpy.lexsort((t["time_start"][m], key))
        m = m[order]
        key = key[order]
        bounds = numpy.flatnonzero(numpy.diff(key)) + 1
        for grp in numpy.split(m, bounds):
            tt = t["time_start"][grp]
            ch = t["channel_peak"][grp]
            adc = t["adc_integral"][grp]
            nb = adc >= neighbor
            boosted = self.ta_boosted[grp]
            dom = self.ta_dom[grp]
            ok = (numpy.abs(tt[:, None] - tt[None, :]) <= dt) & (numpy.abs(ch[:, None] - ch[None, :]) <= dch)
            ok &= nb[None, :]
            ok &= ~boosted[None, :] | (dom[None, :] == dom[:, None])
            numpy.fill_diagonal(ok, True)  # the seed counts once, always
            out[grp] = ok.astype(float) @ adc
        return out

    def sum_selection(self, cut, seed, sums, plane=2):
        t = self.tas
        return self.plane_mask(plane) & (t["adc_integral"] >= seed) & (sums >= cut)

    # measures

    def found(self, sel):
        """Interactions giving at least frac of a selected TA's charge."""
        good = self.link_dominant & sel[self.link_ta] & (self.link_int >= 0)
        found = numpy.zeros(len(self.ints["id"]), dtype=bool)
        found[self.link_int[good]] = True
        return found

    def efficiency(self, sel, found=None):
        if found is None:
            found = self.found(sel)
        nvis = self.visible.sum()
        nfound = (found & self.visible).sum()
        eff = nfound / nvis if nvis else float("nan")
        err = numpy.sqrt(eff * (1 - eff) / nvis) if nvis else float("nan")
        return eff, err

    def coverage(self, sel, block):
        """Covered fraction of the readout by blocks of the selected TAs."""
        t = self.tas
        fr = t["frame"][sel].astype(numpy.int64)
        mid = 0.5 * (t["time_start"][sel] + t["time_end"][sel])
        a = numpy.clip(mid - 0.5 * block, self.lo[fr], self.hi[fr])
        b = numpy.clip(mid + 0.5 * block, self.lo[fr], self.hi[fr])
        shift = fr * 2 * (numpy.max(self.hi - self.lo) + block + 1.0)
        return union_length(a + shift, b + shift) / (self.live * s)

    def weighted_coverage(self, sel, block):
        """Coverage at physical rates: unboosted TAs merge exactly, boosted
        ones are added per boost (sparse approximation)."""
        cov = self.coverage(sel & ~self.ta_boosted, block)
        for b in numpy.unique(self.ta_boost[sel & self.ta_boosted]):
            cov += self.coverage(sel & self.ta_boosted & (self.ta_boost == b), block) / b
        return min(cov, 1.0)

    def full_rate(self, bytes_per_sample):
        """Bytes per second to read all channels of the unit continuously."""
        return self.nchannels / (self.tick / s) * bytes_per_sample

    def metrics(self, sel, block=200 * us, bytes_per_sample=None, relaxed=None):
        """Measures of a selection, at physical rates.

        relaxed, a looser version of sel (eg half the cut), names the
        background sources that could pass: when none passes sel, the 95%
        CL limit is set by the smallest boost among them (sources whose end
        points are far below the cut cannot contribute)."""
        if bytes_per_sample is None:
            bytes_per_sample = self.meta.get("adc_resolution", 14) / 8.0
        w = 1.0 / self.ta_boost
        bkg = sel & ~self.ta_snb & ~self.ta_stress
        sig = sel & self.ta_signal
        rel = sel & self.ta_snb
        r_bkg = w[bkg].sum() / self.live
        r_sig = w[sig].sum() / self.live
        r_rel = w[rel].sum() / self.live
        found = self.found(sel)
        eff, err = self.efficiency(sel, found)
        f_bkg = self.weighted_coverage(bkg, block)
        f_rel = min(self.weighted_coverage(rel, block), 1.0 - f_bkg)
        full = self.full_rate(bytes_per_sample)
        r_all = r_bkg + r_rel
        # a 95% CL limit when no background TA passes, at the smallest
        # boost of the background sources that could pass
        candidates = [i for i, src in enumerate(self.sources)
                      if not src.get("signal", False) and not self.source_stress[i]]
        if relaxed is not None:
            near = set(self.ta_source[relaxed & ~self.ta_snb & ~self.ta_stress & (self.ta_source >= 0)])
            if near:
                candidates = [i for i in candidates if i in near]
            else:  # none near: the MeV-scale sources, those boosted
                candidates = [i for i in candidates if self.source_boost[i] > 1] or candidates
        limit = 3.0 / self.live / min([self.source_boost[i] for i in candidates], default=1.0)
        by_source = {}
        for i, name in enumerate(self.source_names):
            mm = bkg & (self.ta_source == i)
            if mm.any():
                by_source[name] = w[mm].sum() / self.live
        if (bkg & (self.ta_source < 0)).any():
            by_source["noise"] = w[bkg & (self.ta_source < 0)].sum() / self.live
        return dict(
            efficiency=eff, efficiency_err=err,
            rate_signal_hz=r_sig, rate_snb_related_hz=r_rel, rate_background_hz=r_bkg,
            n_background=int(bkg.sum()), n_signal=int(sig.sum()),
            rate_background_limit_hz=(limit if bkg.sum() == 0 else float("nan")),
            background_by_source=by_source,
            purity=r_sig / r_all if r_all > 0 else float("nan"),
            coverage_background=f_bkg, coverage_burst=f_bkg + f_rel,
            data_background_Bps=f_bkg * full, data_burst_Bps=(f_bkg + f_rel) * full,
            data_full_Bps=full,
            pass_probability=self.pass_probability(found),
        )

    def pass_probability(self, found):
        """Per source: probability that one decay (arriving in the readout)
        is found by the selection."""
        out = {}
        for i, name in enumerate(self.source_names):
            m = self.int_source == i
            ndecays = m.sum() * self.window_fraction
            out[name] = (found & m).sum() / ndecays if ndecays > 0 else float("nan")
        return out

    def efficiency_vs(self, sel, column="energy", edges=(0, 10, 15, 20, 25, 30, 40, 60)):
        found = self.found(sel)
        vals = self.ints[column]
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = self.visible & (vals >= lo) & (vals < hi)
            n = m.sum()
            k = (found & m).sum()
            e = k / n if n else float("nan")
            out.append((lo, hi, int(n), e, numpy.sqrt(e * (1 - e) / n) if n else float("nan")))
        return out

    def adc_per_mev(self, plane=2):
        """Collection ADC integral per MeV deposited, from SNB interactions:
        the sum of adc_integral of the TAs each dominates over its edep."""
        t = self.tas
        good = self.link_dominant & (self.link_int >= 0)
        good &= (t["plane"][self.link_ta] == plane) & (self.ints["signal"][numpy.maximum(self.link_int, 0)] == 1)
        adc = numpy.zeros(len(self.ints["id"]))
        numpy.add.at(adc, self.link_int[good], t["adc_integral"][self.link_ta[good]])
        m = self.visible & (adc > 0)
        return float(numpy.median(adc[m] / self.ints["edep"][m]))


def human(bps):
    for unit, scale in (("TB/s", 1e12), ("GB/s", 1e9), ("MB/s", 1e6), ("kB/s", 1e3)):
        if bps >= scale:
            return f"{bps/scale:.3g} {unit}"
    return f"{bps:.3g} B/s"
