---
name: check-dunereco-config
description: Audit whether an SPNG signal-processing configuration stays faithful to the official dunereco OmnibusSigProc (OSP) config, and propose compliant changes. Use when asked to check/verify SPNG config faithfulness, compare SPNG vs official/dunereco/OSP signal processing, or when the SPNG cfg under toolkit/spng/cfg/spng may have drifted from reference/dunereco.
---

# Checking SPNG config faithfulness to official dunereco OSP

The SPNG configuration under `toolkit/spng/cfg/spng` was derived from the
"official" dunereco configuration that drives WCT's original OmnibusSigProc
(OSP). SPNG evaluations must stay faithful to that official config. The
authoritative source is the dunereco reference tree, by default
`reference/dunereco/dunereco/DUNEWireCell/{pdhd,pdsp,protodunevd}/`.

The tool `toolkit/spng/test/check-dunereco-config/check-dunereco-config` drives
both halves of the check. Always run it from a shell where the WCT build is
available (it uses `install/bin/wcsonnet`). It evaluates each Jsonnet config to
JSON with the required TLAs/ext-vars, so you do not have to construct those.

## Comparisons (run `list` first)

Pick a comparison with `-c` (or `--detector`). Run `./check-dunereco-config
list` to see them. The important subtlety is **PDHD's bad first APA**: the real
detector's first APA is partly broken with a non-standard field response, and
the official PDHD config bakes in bad-APA tunings (`plane2layer: [0,2,1]` plane
swap, per-APA `Wiener_tight_*_APA1` filters). SPNG idealises PDHD as all-good.
So `-c pdhd` (SPNG PDHD vs official PDHD) is **confounded** by those quirks.

- `-c pdhd-x-pdsp` is the **cross comparison**: SPNG PDHD vs official **PDSP**,
  the nominal all-good ProtoDUNE-SP design (same layout, no bad-APA tunings). It
  is the faithful reference for an idealised SPNG PDHD. Prefer it for judging
  whether SPNG PDHD is correct.
- In a cross comparison some deltas are genuine **detector** differences, not
  SPNG infidelities (e.g. `ADC_mV` differs ~4x: PDSP is 12-bit, PDHD 14-bit).
  Call these out separately from real infidelities.
- `-c pdvd` compares SPNG PDVD vs official protodunevd. PDVD has two drift
  volumes (tpcid 0..3 = bottom, 4..7 = top) with different electronics/ADC; the
  official fixture uses a bottom anode, so the comparison uses a bottom tpcid.
- `-c pdsp` is currently unavailable (no SPNG PDSP config). The tool reports it
  clearly; use `pdhd-x-pdsp` for the nominal check.

## Two-part method

**Part A — mechanistic diff (run first).** A filtered, like-for-like JSON diff
of the evaluated OSP configuration: the official `OmnibusSigProc` + `HfFilter`
/`LfFilter` nodes vs the SPNG "mirror OSP" config (driven through
`spng/adc-to-osp.jsonnet`). Superfluous differences (instance names, per-anode
tag suffixes, node-handle instance names) are normalised away.

```
toolkit/spng/test/check-dunereco-config/check-dunereco-config diff -c pdhd
toolkit/spng/test/check-dunereco-config/check-dunereco-config diff -c pdhd-x-pdsp
```

Each finding is printed with a concrete `fix:` line. `DRIFT`/`MISSING` are
genuine parameter/filter deltas; `EXTRA` are SPNG-only keys. Exit code is
non-zero when drift is found. Use `--json` for machine-readable output and
`--all-diffs` to also show deltas the registry marks as intentional.

**Part B — semantic comparison (LLM; that's you).** The SPNG-native "tensor data
model" (TDM) config expresses the same physics in a different idiom
(`SPNGFilterKernel`/`SPNGResponseKernel`/`SPNGDeconKernel`/`SPNGKernelConvolve`),
so a syntactic diff cannot map it. Emit the self-contained prompt:

```
toolkit/spng/test/check-dunereco-config/check-dunereco-config emit-prompt -c pdhd-x-pdsp -o /tmp/cdc-pdhd.md
```

Read that file. It contains the official OSP parameters + filters, the
SPNG-native TDM node extracts, and a mapping rubric. Work through the rubric to
decide, for each official filter and OSP parameter, whether the SPNG-native
config reproduces it. Then produce a compliance report.

## What to deliver

Produce a single report with two sections:

1. **Mechanistic findings (Part A):** relay the `diff` findings, and for each,
   classify it as **infidelity** (must fix) or **intentional SPNG adaptation**
   (explain why). Use judgement: e.g. an APA0-specific special-case that the
   SPNG mirror omits is a real infidelity for APA0; a debug/dump tag that only
   makes sense in one runtime is not. When unsure, flag it rather than dismiss.

2. **Semantic findings (Part B):** the mapped comparison from the emitted
   prompt — official element, SPNG-native equivalent, match/mismatch, severity,
   and the minimal Jsonnet edit (file + change) to fix each mismatch. Point at
   real files under `toolkit/spng/cfg/spng` (e.g. `detconfs/<det>.jsonnet` for
   filter magic numbers, `detconfs/<det>/sp-filters.jsonnet` and
   `detconfs/<det>/sp.jsonnet` for the mirror OSP config, `decon.jsonnet` for
   response-kernel scaling).

Offer the changes; do not apply them unless the user asks. If asked to fix, make
the minimal edit and re-run `diff` to confirm the finding clears.

## Adding a detector / comparison

The registry at the top of the script has three parts: `OFFICIAL[det]` (fixture
template + anode + extvars), `SPNG[det]` (the adc-to-osp/adc-to-spng jobs, or
`None` if SPNG lacks that detector), and `COMPARISONS[name]` (a `(spng,
official)` pairing). To add an official detector, add a
`fixtures/official-<det>.jsonnet` (copy the pdhd one, adjust the
`import '<det>/sp.jsonnet'` line, params/tools and the make_sigproc call) and an
`OFFICIAL` entry. To add a SPNG detector, add a `SPNG` entry
(`_spng_jobs("<det>")`). Then add `COMPARISONS` entries. Known gaps: SPNG has no
`pdsp` config and its `pdvd` config currently fails to evaluate. See README.org.

## Gotchas

- The official config needs `-V elecGain=<mV/fC>` (recorded per-case); the SPNG
  jobs need `-A input=...` (a dummy filename is fine — nothing is read at
  evaluation time). The script supplies these.
- The reference tree is imported by absolute path (not via `WIRECELL_PATH`) to
  avoid shadowing shared Jsonnet files; `--reference DIR` overrides its location.
- Temp fixtures go to `$TMPDIR`; set it to a writable dir if `/tmp` is small.
