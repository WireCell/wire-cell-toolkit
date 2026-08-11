#!/usr/bin/env python3
"""Make ROI-wise efficiency/purity plots from a .pt file holding 'y' and 'labels'.

    python plot_roi_metrics.py xvu-test_line_uplane_t1-87-t2-87.pt
    python plot_roi_metrics.py <file>.pt --show          # interactive window too
    python plot_roi_metrics.py <file>.pt --show --no-save

Produces <stem>_roi_metrics.png with four panels:
  1. efficiency & purity vs threshold
  2. purity vs efficiency (the working-point curve, threshold annotated)
  3. efficiency vs true ROI width, purity vs reco ROI width (at --threshold)
  4. ROI width distributions, true vs reco (at --threshold)
"""

import argparse
from pathlib import Path

import matplotlib
import torch

from roi_metrics import roi_metrics, roi_table, threshold_scan


def select_gui_backend():
    """Activate the first usable interactive backend; False if there is none."""
    import importlib
    import os
    import sys

    if sys.platform not in ("darwin", "win32") and not (
            os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        return False
    # matplotlib.use() does not check that the toolkit is importable (it only
    # validates the name until pyplot loads), so import the backend ourselves.
    for name in ("TkAgg", "QtAgg", "GTK4Agg", "WXAgg", "MacOSX"):
        try:
            importlib.import_module(f"matplotlib.backends.backend_{name.lower()}")
        except Exception:
            continue
        matplotlib.use(name)
        return True
    return False


def binned_rate(width, matched, edges):
    """Matched fraction in each width bin, plus binomial errors and bin counts."""
    idx = torch.bucketize(width.to(torch.float64), edges[1:-1], right=True)
    nbin = len(edges) - 1
    n = torch.bincount(idx, minlength=nbin).to(torch.float64)
    k = torch.bincount(idx[matched], minlength=nbin).to(torch.float64)
    p = torch.where(n > 0, k / n, torch.full_like(n, float("nan")))
    err = torch.where(n > 0, ((p * (1 - p)).clamp(min=0) / n).sqrt(), torch.zeros_like(n))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, p, err, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", type=Path, help=".pt file with 'y' and 'labels' tensors")
    ap.add_argument("-t", "--threshold", type=float, default=0.5,
                    help="working-point threshold on y (default: 0.5)")
    ap.add_argument("-o", "--out", type=Path, default=None, help="output image path")
    ap.add_argument("--nscan", type=int, default=41, help="points in the threshold scan")
    ap.add_argument("--show", action="store_true",
                    help="open an interactive window (needs a display; over ssh use -X/-Y)")
    ap.add_argument("--no-save", dest="save", action="store_false",
                    help="skip writing the image (only useful with --show)")
    args = ap.parse_args()

    # Pick the backend before pyplot is imported: an interactive one only when
    # asked for, otherwise headless Agg.
    if args.show:
        if not select_gui_backend():
            print("warning: no interactive backend available "
                  "(no display or no GUI toolkit) -- falling back to a saved image")
            args.show, args.save = False, True
    if not args.show:
        matplotlib.use("Agg")
    global plt
    import matplotlib.pyplot as plt

    d = torch.load(args.path, map_location="cpu")
    y, labels = d["y"], d["labels"]
    print(f"{args.path.name}: y {tuple(y.shape)}  labels {tuple(labels.shape)}")

    wp = roi_metrics(y, labels, threshold=args.threshold)
    print(f"threshold {args.threshold}")
    print(f"  efficiency {float(wp['efficiency']):.4f} "
          f"({int(wp['n_true_matched'])}/{int(wp['n_true'])} true ROIs)")
    print(f"  purity     {float(wp['purity']):.4f} "
          f"({int(wp['n_reco_matched'])}/{int(wp['n_reco'])} reco ROIs)")

    # Scan interior thresholds only; 0 and 1 are degenerate for probabilities.
    thr = torch.linspace(0, 1, args.nscan)[1:-1]
    scan = threshold_scan(y, labels, thr)

    true_tab = roi_table(labels == 1, y > args.threshold)
    reco_tab = roi_table(y > args.threshold, labels == 1)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

    ax = axes[0, 0]
    ax.plot(scan["threshold"], scan["efficiency"], "-o", ms=3, label="efficiency")
    ax.plot(scan["threshold"], scan["purity"], "-s", ms=3, label="purity")
    ax.axvline(args.threshold, color="0.6", ls="--", lw=1)
    ax.set_xlabel("threshold on y")
    ax.set_ylabel("ROI-wise rate")
    ax.set_ylim(0, 1.02)
    ax.set_title("Efficiency and purity vs threshold")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(scan["efficiency"], scan["purity"], "-", lw=1, color="0.5", zorder=1)
    sc = ax.scatter(scan["efficiency"], scan["purity"], c=scan["threshold"],
                    cmap="viridis", s=18, zorder=2)
    ax.scatter([float(wp["efficiency"])], [float(wp["purity"])], marker="*", s=200,
               color="crimson", zorder=3, label=f"t = {args.threshold}")
    fig.colorbar(sc, ax=ax, label="threshold")
    ax.set_xlabel("efficiency")
    ax.set_ylabel("purity")
    ax.set_title("Working-point curve")
    ax.legend(loc="lower left")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    widths = torch.cat([true_tab["length"], reco_tab["length"]]).to(torch.float64)
    if widths.numel():
        # log-spaced: narrow ROIs are where the interesting inefficiency and
        # impurity live, so they deserve their own bins.
        hi = max(float(widths.max()), 2.0)
        edges = torch.logspace(0, torch.log10(torch.tensor(hi)).item(), 13,
                               dtype=torch.float64)
        edges[0] = 0.5
        for tab, name, style in ((true_tab, "efficiency (vs true width)", "-o"),
                                 (reco_tab, "purity (vs reco width)", "-s")):
            c, p, e, n = binned_rate(tab["length"], tab["matched"], edges)
            keep = n > 0
            ax.errorbar(c[keep], p[keep], yerr=e[keep], fmt=style, ms=4,
                        capsize=2, label=name)
    ax.set_xlabel("ROI width [time bins]")
    ax.set_ylabel("matched fraction")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Rates vs ROI width (t = {args.threshold})")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    if widths.numel():
        bins = torch.linspace(0.5, float(widths.max()) + 0.5,
                              min(40, int(widths.max()) + 1)).tolist()
        ax.hist(true_tab["length"].tolist(), bins=bins, histtype="step", lw=1.6,
                label=f"true ({len(true_tab['length'])})")
        ax.hist(reco_tab["length"].tolist(), bins=bins, histtype="step", lw=1.6,
                label=f"reco ({len(reco_tab['length'])})")
    ax.set_xlabel("ROI width [time bins]")
    ax.set_ylabel("ROIs")
    ax.set_yscale("log")
    ax.set_title("ROI width distributions")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.suptitle(f"ROI metrics — {args.path.name}", fontsize=12)
    fig.tight_layout()

    if args.save:
        out = args.out or args.path.with_name(args.path.stem + "_roi_metrics.png")
        fig.savefig(out, dpi=140)
        print(f"wrote {out}")
    if args.show:
        print(f"showing interactively ({matplotlib.get_backend()}); close the window to exit")
        plt.show()


if __name__ == "__main__":
    main()
