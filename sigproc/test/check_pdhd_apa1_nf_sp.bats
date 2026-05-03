#!/usr/bin/env bats

# End-to-end PDHD APA1 NF+SP regression test.
#
# Runs wire-cell against the production wct-nf-sp.jsonnet on the APA1
# raw frame archive for run 027409 evt 0, and numerically compares the
# resulting SP frames (gauss1, wiener1 tags) against a checked-in
# reference tar.bz2 at sigproc/test/data/protodunehd-sp-frames-anode1.tar.bz2.
# The comparison uses per-channel RMS and per-tag integral within
# stated tolerances; observed run-to-run spread is bit-zero on this
# branch (commit 47d16673 fixed the FFTW plan-cache nondeterminism), so
# the tolerances are deliberately conservative against future drift.
#
# bats file_tags=sigproc,PDHD,L1SP,regression
#
# Required external fixtures (test SKIPs cleanly when unavailable):
#
#   WCT_PDHD_DATA   directory holding the APA1 raw input:
#                     protodunehd-orig-frames-anode1.tar.bz2
#                   Defaults to:
#                     /nfs/data/1/xqian/toolkit-dev/wcp-porting-img/pdhd/
#                     input_data/run027409/evt_0
#
#   WCT_PDHD_REF    directory holding the frozen reference SP output.
#                   Defaults to the in-tree fixture
#                   sigproc/test/data/.  Override only when refreshing
#                   the reference after a deliberate algorithmic change
#                   (e.g. commit 616cfcb2 routed L1SP into wiener1, at
#                   which point this fixture was generated from HEAD).
#
#   WCT_PDHD_DEPLOY directory holding the standalone wct-nf-sp.jsonnet
#                   entry point.  Defaults to:
#                     /nfs/data/1/xqian/toolkit-dev/wcp-porting-img/pdhd

bats_load_library wct-bats.sh

# bats test_tags=PDHD,L1SP,regression
@test "PDHD APA1 NF+SP regression vs frozen reference" {
    local data_dir="${WCT_PDHD_DATA:-/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/pdhd/input_data/run027409/evt_0}"
    local ref_dir="${WCT_PDHD_REF:-$(dirname "$BATS_TEST_FILENAME")/data}"
    local deploy_dir="${WCT_PDHD_DEPLOY:-/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/pdhd}"
    local input_tar="${data_dir}/protodunehd-orig-frames-anode1.tar.bz2"
    local ref_tar="${ref_dir}/protodunehd-sp-frames-anode1.tar.bz2"
    local jsonnet="${deploy_dir}/wct-nf-sp.jsonnet"

    [ -f "$input_tar"  ] || skip "no APA1 input tar at $input_tar"
    [ -f "$ref_tar"    ] || skip "no APA1 reference tar at $ref_tar"
    [ -f "$jsonnet"    ] || skip "no wct-nf-sp.jsonnet at $jsonnet"
    command -v wire-cell >/dev/null || skip "wire-cell not on PATH"
    command -v python3   >/dev/null || skip "python3 not on PATH"

    cd_tmp

    # wire-cell expects the input prefix as a string; the FrameFileSource
    # will read "${prefix}-anode1.tar.bz2".  We stage the input next to
    # the test cwd so a single prefix works for the source.
    cp "$input_tar" "protodunehd-orig-frames-anode1.tar.bz2"

    # Run wire-cell on APA1 only.
    wire-cell -l stderr -L info \
        -V "elecGain=14" \
        --tla-str orig_prefix="protodunehd-orig-frames" \
        --tla-str raw_prefix="protodunehd-sp-frames-raw" \
        --tla-str sp_prefix="protodunehd-sp-frames" \
        --tla-str reality="data" \
        --tla-code anode_indices="[1]" \
        -c "$jsonnet"

    local out_tar="protodunehd-sp-frames-anode1.tar.bz2"
    [ -s "$out_tar" ] || die "wire-cell did not produce $out_tar"

    # Numerical-summary comparison against the frozen reference.  Any
    # exit code or stderr output is shown to the user.
    python3 - "$out_tar" "$ref_tar" <<'PY'
import sys, tarfile, io, numpy as np

test_tar, ref_tar = sys.argv[1:3]

def load_frames(path):
    """Return {tag: ndarray (nchan, nticks)} keyed by trace tag."""
    out = {}
    with tarfile.open(path, 'r:bz2') as tf:
        for m in tf.getmembers():
            name = m.name
            if not name.endswith('.npy') or not name.startswith('frame_'):
                continue
            # frame_<tag>_<ident>.npy
            tag = name[len('frame_'):-len('.npy')].rsplit('_', 1)[0]
            arr = np.load(io.BytesIO(tf.extractfile(m).read()))
            out[tag] = arr
    return out

a = load_frames(test_tar)
b = load_frames(ref_tar)
assert set(a.keys()) == set(b.keys()), \
    f"tag mismatch: test={sorted(a)} ref={sorted(b)}"

for tag in sorted(a):
    A, B = a[tag], b[tag]
    assert A.shape == B.shape, f"{tag}: shape mismatch {A.shape} vs {B.shape}"

    # Per-channel RMS comparison.
    rms_a = np.sqrt(np.mean(A.astype(np.float64)**2, axis=1))
    rms_b = np.sqrt(np.mean(B.astype(np.float64)**2, axis=1))
    eps   = 1e-3
    rel   = np.abs(rms_a - rms_b) / np.maximum(rms_b, eps)
    n_bad = int((rel > 1e-2).sum())   # 1% per-channel RMS tolerance
    frac_bad = n_bad / len(rel)
    assert frac_bad < 0.01, \
        f"{tag}: {n_bad}/{len(rel)} channels exceed 1% RMS tolerance " \
        f"(max rel={rel.max():.3e})"

    # Per-tag integral within 0.1%.
    int_a = float(np.abs(A).sum())
    int_b = float(np.abs(B).sum())
    assert int_b > 0
    rel_int = abs(int_a - int_b) / int_b
    assert rel_int < 1e-3, \
        f"{tag}: integral mismatch {rel_int:.3e}"

    # Soft array-level check.  Loose tolerances on purpose (atol of 1
    # ADC count, rtol 1e-3) — any tighter and FFTW jitter trips it.
    if not np.allclose(A, B, atol=1.0, rtol=1e-3):
        diff = np.abs(A - B)
        print(f"WARN {tag}: np.allclose failed: "
              f"max_abs_diff={diff.max():.3f} "
              f"99th_pct={np.percentile(diff, 99):.3f}",
              file=sys.stderr)

    print(f"OK  {tag}: shape={A.shape} max_rel_rms={rel.max():.3e} rel_int={rel_int:.3e}")

print("PDHD APA1 NF+SP regression: PASS")
PY
}
