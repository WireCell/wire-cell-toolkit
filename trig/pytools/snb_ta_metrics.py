#!/usr/bin/env python3
'''
Efficiency, purity and readout data rate of trigger activities (TAs) as
seeds for supernova neutrino burst (SNB) readout, from the truth tables of
trig/cfg/snb-ta-study.jsonnet.  See snbta.py for the definitions.

Rates are per simulated unit (HD APA or VD CRU) and at the physical SNB rate
(the boost is divided out).  "module" scales by the module's active mass
over the unit's (sparse approximation).
'''

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import snbta  # noqa: E402
from snbta import human, us, cm  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="TATruth .npz files")
    p.add_argument("--plane", default="2", help="plane index of seed TAs, or 'all' (default 2, collection)")
    p.add_argument("--frac", type=float, default=0.5, help="charge fraction for signal TA and for found")
    p.add_argument("--margin-us", type=float, default=50.0, help="visible: margin from readout ends")
    p.add_argument("--block-us", type=float, default=200.0, help="readout block length per TA")
    p.add_argument("--module-kt", type=float, default=14.0, help="active argon mass of a module")
    p.add_argument("--charge-cuts", default="0,1000,2000,4000,6000,8000,10000,12000,16000,20000",
                   help="method 'charge': TA adc_integral cuts to scan")
    p.add_argument("--sum-cuts", default="10000,20000,30000,40000,50000,60000",
                   help="method 'sum': neighborhood sum cuts to scan")
    p.add_argument("--sum-seed", type=float, default=3000, help="method 'sum': seed TA adc_integral")
    p.add_argument("--sum-neighbor", type=float, default=1000, help="method 'sum': neighbor TA adc_integral")
    p.add_argument("--sum-radius-cm", type=float, default=30.0, help="method 'sum': neighborhood radius")
    p.add_argument("--json", default=None, help="write all metrics to this JSON file")
    args = p.parse_args(argv)

    st = snbta.Study(args.files, frac=args.frac, margin_us=args.margin_us)
    plane = None if args.plane == "all" else int(args.plane)
    module = args.module_kt * 1e6 / st.mass_kg
    m = st.meta
    print(f"{m.get('geometry', st.detector)}: {st.mass_kg/1000:.1f} t, {st.nchannels} channels, "
          f"{st.nframes} frames, live {st.live:.2f} s")
    print(f"TA {m['ta_algorithm']}, thresholds {m['thresholds']} ADC, noise {m['noise']}; "
          f"SNB {st.snb_rate_hz:g} Hz simulated, {st.snb_physical_rate_hz:.3g} Hz physical, boost {st.boost:.0f}")
    print(f"visible SNB interactions {st.visible.sum()}; collection ADC per MeV {st.adc_per_mev():.0f}")
    print()

    hdr = (f"{'method':>8} {'cut':>7} {'eff':>12} {'sig/s':>7} {'bkg/s':>9} {'purity':>8} "
           f"{'data bkg':>10} {'burst':>10} {'module bkg':>11}")
    results = []

    def row(method, cut, sel):
        r = st.metrics(sel, block=args.block_us * us)
        r.update(method=method, cut=cut)
        results.append(r)
        bkg = f"{r['rate_background_hz']:9.3g}" if r["n_background"] else f"<{r['rate_background_limit_hz']:8.2g}"
        print(f"{method:>8} {cut:7.0f} {r['efficiency']:6.3f}±{r['efficiency_err']:.3f} {r['rate_signal_hz']:7.3g} "
              f"{bkg} {r['purity']:8.3f} {human(r['data_background_Bps']):>10} {human(r['data_burst_Bps']):>10} "
              f"{human(r['data_background_Bps'] * module):>11}")

    print(hdr)
    for cut in [float(x) for x in args.charge_cuts.split(",")]:
        row("charge", cut, st.charge_selection(cut, plane=plane))
    sums = st.neighborhood_sum(plane=plane, radius=args.sum_radius_cm * cm, neighbor=args.sum_neighbor)
    for cut in [float(x) for x in args.sum_cuts.split(",")]:
        row("sum", cut, st.sum_selection(cut, args.sum_seed, sums, plane=plane))
    print()
    print(f"sig/s, bkg/s: physical rates of signal and background TAs per unit; full readout "
          f"{human(results[0]['data_full_Bps'])} per unit; module = {module:.0f} units.")
    print("per-decay pass probability of the background sources, for each row:")
    names = [n for n in st.source_names if n != "snb"]
    print(f"{'method':>8} {'cut':>7} " + " ".join(f"{n:>9}" for n in names))
    for r in results:
        print(f"{r['method']:>8} {r['cut']:7.0f} " + " ".join(f"{r['pass_probability'][n]:9.2e}" for n in names))

    if args.json:
        with open(args.json, "w") as fp:
            json.dump(dict(meta=m, nframes=st.nframes, live_s=st.live, boost=st.boost, results=results,
                           args=vars(args)), fp, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
