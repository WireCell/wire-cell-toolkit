#!/usr/bin/env python3
"""Self-test for check-dunereco-config.

Two layers:

  * Pure-logic tests of the normalisation/diff engine -- no build, no wcsonnet
    required.  These always run.

  * An optional live smoke test that actually evaluates the pdhd configs and
    asserts the diff engine returns well-formed structure.  Skipped (not failed)
    when wcsonnet / the reference tree are unavailable.

Run:  python3 selftest.py
Exit: 0 on success, 1 on failure.
"""

import importlib.util
import os
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_module():
    # The tool has no .py extension; load it by path with an explicit loader.
    path = HERE / "check-dunereco-config"
    loader = SourceFileLoader("cdc", str(path))
    spec = importlib.util.spec_from_loader("cdc", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


FAILS = []


def check(cond, msg):
    if cond:
        print("  ok  : %s" % msg)
    else:
        print("  FAIL: %s" % msg)
        FAILS.append(msg)


def test_pure(cdc):
    print("[pure] normalisation")
    # Node-handle keys: compare type only, keep "" meaningful.
    check(cdc.norm_handle("FieldResponse:pdhd_fr") == "FieldResponse",
          "node handle -> type token")
    check(cdc.norm_handle("") == "", "empty handle stays empty")

    # Tag suffix: strip trailing per-anode digits.
    check(cdc.strip_tag_suffix("wiener3") == "wiener", "tag suffix stripped")
    check(cdc.strip_tag_suffix("") == "", "empty tag stays empty")

    d = cdc.normalize_osp_data({
        "anode": "AnodePlane:apa7",
        "per_chan_resp": "",
        "wiener_tag": "wiener2",
        "lroi_rebin": 6,
    })
    check(d["anode"] == "AnodePlane", "normalize_osp_data handle")
    check(d["per_chan_resp"] == "", "normalize_osp_data empty handle kept")
    check(d["wiener_tag"] == "wiener", "normalize_osp_data tag stripped")
    check(d["lroi_rebin"] == 6, "normalize_osp_data plain value kept")

    print("[pure] value equality (float tolerance)")
    check(cdc.values_equal(1.0, 1.0 + 1e-15), "close floats equal")
    check(not cdc.values_equal(2e-6, 3e-6), "tau 2e-6 != 3e-6 (the drift case)")
    check(cdc.values_equal([1, 2, 3], [1, 2, 3]), "equal lists")
    check(not cdc.values_equal([1, 2], [1, 2, 3]), "unequal-length lists differ")
    check(not cdc.values_equal(True, 1), "bool not conflated with int")

    print("[pure] diff_data")
    off = {"a": 1, "b": 2, "only_off": 9}
    spng = {"a": 1, "b": 5, "only_spng": 8}
    dd = cdc.diff_data(off, spng)
    check(dd["only_official"] == ["only_off"], "only_official detected")
    check(dd["only_spng"] == ["only_spng"], "only_spng detected")
    check(dd["changed"] == [("b", 2, 5)], "changed value detected")

    print("[pure] comparison registry")
    check("pdhd-x-pdsp" in cdc.COMPARISONS, "cross comparison registered")
    check(cdc.COMPARISONS["pdhd-x-pdsp"]["spng"] == "pdhd"
          and cdc.COMPARISONS["pdhd-x-pdsp"]["official"] == "pdsp",
          "pdhd-x-pdsp pairs SPNG pdhd with official pdsp")
    check(cdc.SPNG.get("pdsp") is None, "SPNG pdsp correctly absent")
    # Unavailable SPNG side must raise a clear EvalError (no eval needed).
    try:
        cdc.resolve_comparison("pdsp")
        check(False, "resolve_comparison('pdsp') should raise")
    except cdc.EvalError as exc:
        check("no 'pdsp' configuration" in str(exc),
              "unavailable comparison raises explanatory EvalError")


def test_live(cdc):
    print("[live] smoke test (needs wcsonnet + reference tree)")
    wcsonnet = cdc.default_wcsonnet()
    reference = Path(os.environ.get("CDC_REFERENCE", cdc.default_reference()))
    have_wc = wcsonnet != "wcsonnet" and Path(wcsonnet).exists()
    if not (have_wc and reference.exists()):
        print("  skip: wcsonnet=%s exists=%s  reference exists=%s"
              % (wcsonnet, have_wc, reference.exists()))
        return
    try:
        result = cdc.compare("pdhd", wcsonnet, reference)
    except cdc.EvalError as exc:
        print("  skip: evaluation unavailable in this environment:\n    %s"
              % str(exc).splitlines()[0])
        return
    check("osp" in result and "changed" in result["osp"], "diff produced osp block")
    check(isinstance(result["filters_changed"], list), "filters_changed is a list")
    # The engine must faithfully surface a filter-magic-number difference when
    # one exists.  We do not assert a *specific* current drift (that would
    # invert the moment it is fixed); we assert the machinery ran end-to-end and
    # the ADC_mV / filter comparison paths are exercised.
    keys = set(result["osp_official"]) & set(result["osp_spng"])
    check("ADC_mV" in keys, "ADC_mV present on both sides (comparison exercised)")


def main():
    cdc = load_module()
    test_pure(cdc)
    test_live(cdc)
    print()
    if FAILS:
        print("SELFTEST FAILED: %d check(s)" % len(FAILS))
        return 1
    print("SELFTEST PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
