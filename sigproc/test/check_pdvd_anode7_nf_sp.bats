#!/usr/bin/env bats

# End-to-end PDVD top-drift anode 7 NF+SP regression test.
#
# Runs wire-cell against a pre-rendered frozen config on the anode 7
# raw frame archive for run 039324 evt 0, and compares the resulting SP
# frames (gauss7, wiener7 tags) bit-exactly against the checked-in
# reference tar at sigproc/test/data/protodunevd-sp-frames-anode7.tar.bz2.
#
# Differences from the PDHD APA1 test:
#   - reality='data': on top anodes (ident >= 4) this flag is a no-op;
#     the 512->500 ns Resampler is gated on n<4 in pdvd/wct-nf-sp.jsonnet.
#   - l1sp_pd_mode='process': LASSO writeback active (production top path;
#     the former process->dump auto-downgrade for ident >= 4 was removed).
#   - No elecGain override — protodunevd/params.jsonnet carries its own gains.
#
# Determinism: SP output (wiener7, gauss7) is bit-identical run-to-run on
# this branch after two fixes:
#
#   1. Use-after-free + iterator UB + double-free in BreakROIs/BreakROI1
#      iterating rois_*_loose while BreakROI erases and deletes from it.
#      Snapshot-then-iterate pattern; BreakROI1's second-pass loop disabled
#      because it operated on the freed originals.  See ROI_refinement.cxx.
#
#   2. Uninitialized FFT-padding rows in OmnibusSigProc::decon_2D_looseROI.
#      That function only populated c_data_afterfilter rows iterated by
#      m_channel_range[plane] (OSP wires 0..m_nwires-1) and left rows
#      m_nwires..m_fft_nwires-1 uninitialized.  inv_c2r turned that into
#      heap garbage, which the m_pad_nwires-offset .block() extract pulled
#      into the LAST m_pad_nwires rows of m_r_data — for PDVD V plane
#      (m_pad_nwires[1]=11), that was OSP wires 465-475 = WCT idents
#      1228-1238 = frame rows 752-762.  The corrupted m_r_data fed
#      find_ROI_loose, producing ASLR-dependent ROI bounds that survived
#      all downstream ROI refinement and ended up written back into
#      m_r_data via apply_roi.  Fix mirrors decon_2D_tightROI /
#      decon_2D_ROI_refine / decon_2D_charge: initialize all rows with
#      the default filter, then override per-channel rows that have
#      bad/lf_noisy neighbors.  See OmnibusSigProc.cxx::decon_2D_looseROI.
#
# Comparison is bit-exact (np.array_equal on every frame_* npy).  If the
# reference fixture and the test output diverge after a deliberate
# algorithmic change, refresh the reference by re-running this test's
# wire-cell command and replacing the in-tree tar.bz2.
#
# bats file_tags=sigproc,PDVD,L1SP,regression
#
# Required external fixtures (test SKIPs cleanly when unavailable):
#
#   WCT_PDVD_DATA   directory holding the anode 7 raw input:
#                     protodune-orig-frames-anode7.tar.bz2
#                   Defaults to:
#                     /nfs/data/1/xqian/toolkit-dev/toolkit/pdvd/
#                     input_data/run039324/evt_0
#
#   WCT_PDVD_REF    directory holding the frozen reference SP output.
#                   Defaults to the in-tree fixture
#                   sigproc/test/data/.  Override only when refreshing
#                   the reference after a deliberate algorithmic change.
#
# The wire-cell configuration is taken from the in-tree pre-rendered
# snapshot sigproc/test/data/pdvd-anode7-nf-sp-cfg.json (frozen at
# jsonnet render time).  To refresh after a deliberate config change,
# re-run:
#   jsonnet -J cfg \
#     --tla-str orig_prefix=protodune-orig-frames \
#     --tla-str raw_prefix=protodune-sp-frames-raw \
#     --tla-str sp_prefix=protodune-sp-frames \
#     --tla-str reality=data \
#     --tla-code 'anode_indices=[7]' \
#     --tla-str l1sp_pd_mode=process \
#     pdvd/wct-nf-sp.jsonnet > sigproc/test/data/pdvd-anode7-nf-sp-cfg.json
# and commit the new JSON together with the refreshed reference tar.

bats_load_library wct-bats.sh

# bats test_tags=PDVD,L1SP,regression
@test "PDVD anode 7 NF+SP regression vs frozen reference" {
    local data_dir="${WCT_PDVD_DATA:-/nfs/data/1/xqian/toolkit-dev/toolkit/pdvd/input_data/run039324/evt_0}"
    local ref_dir="${WCT_PDVD_REF:-$(dirname "$BATS_TEST_FILENAME")/data}"
    local frozen_cfg="$(dirname "$BATS_TEST_FILENAME")/data/pdvd-anode7-nf-sp-cfg.json"
    local input_tar="${data_dir}/protodune-orig-frames-anode7.tar.bz2"
    local ref_tar="${ref_dir}/protodunevd-sp-frames-anode7.tar.bz2"

    [ -f "$input_tar"  ] || skip "no anode 7 input tar at $input_tar"
    [ -f "$ref_tar"    ] || skip "no anode 7 reference tar at $ref_tar"
    [ -f "$frozen_cfg" ] || skip "no frozen cfg at $frozen_cfg"
    command -v wire-cell >/dev/null || skip "wire-cell not on PATH"
    command -v python3   >/dev/null || skip "python3 not on PATH"

    cd_tmp

    cp "$input_tar" "protodune-orig-frames-anode7.tar.bz2"

    # Run wire-cell using the pre-rendered frozen config (no live jsonnet).
    wire-cell -l stderr -L info -c "$frozen_cfg"

    local out_tar="protodune-sp-frames-anode7.tar.bz2"
    [ -s "$out_tar" ] || die "wire-cell did not produce $out_tar"

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
    if not np.array_equal(A, B):
        # Helpful failure detail: which rows differ and by how much.
        eq = np.all(A == B, axis=1)
        diff_rows = np.where(~eq)[0].tolist()
        d = np.abs(A[diff_rows] - B[diff_rows]) if diff_rows else np.array([0.0])
        raise AssertionError(
            f"{tag}: not bit-exact ({len(diff_rows)} rows differ, "
            f"max|d|={float(d.max()):.4g}, indices[:10]={diff_rows[:10]})")
    print(f"OK  {tag}: shape={A.shape} bit-exact vs reference")

print("PDVD anode 7 NF+SP regression: PASS")
PY
}
