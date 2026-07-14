// Kernel-gated tests for L1SPFilterPD that exercise the lazy
// init_resp() path against the production kernel
// (pdhd_l1sp_kernels.json.bz2, found via WIRECELL_PATH).
//
// Two scopes:
//   1. Kernel-load smoke: configure with a real kernel, run operator()
//      on a minimal frame, assert no exception and well-formed output.
//   2. Dump-mode NPZ schema pin: same input but with dump_mode=true,
//      assert the documented set of frame-level + per-ROI columns is
//      written.  Per-ROI columns may be empty (no triggered ROIs on a
//      minimal flat input) — only the schema is pinned, not values.
//
// File is named check_*.cxx — currently doctest discovery globs
// `doctest*.cxx` so this file does NOT auto-build into wcdoctest-sigproc.
// To enable, rename to doctest_*.cxx (kernel availability is checked at
// run-time via SUBCASE skip, so it is safe to mix into the default
// doctest binary).
//
// Skipping behavior: if WIRECELL_PATH is unset or the kernel file does
// not resolve, every TEST_CASE prints a single SKIP line and exits
// without failing.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/cnpy.h"

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include <cstdlib>
#include <filesystem>
#include <set>
#include <string>

using namespace WireCell;
using namespace WireCell::Aux;

namespace {

const std::string KERNEL = "pdhd_l1sp_kernels.json.bz2";

bool kernel_available()
{
    try {
        auto path = Persist::resolve(KERNEL);
        return !path.empty();
    } catch (...) {
        return false;
    }
}

void load_plugins()
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");
    pm.add("WireCellSigProc");
}

// Build a minimal frame with `raw` and `gauss` tags.  4 channels,
// 1000 ticks of zeros — enough to drive operator() through the full
// per-channel loop without producing any triggered ROIs.
IFrame::pointer make_minimal_frame()
{
    const int nchan = 4;
    const int nticks = 1000;
    ITrace::vector adc_traces, gauss_traces;
    IFrame::trace_list_t raw_idx, gauss_idx;
    for (int ch = 0; ch < nchan; ++ch) {
        Waveform::realseq_t z(nticks, 0.0f);
        adc_traces.push_back(std::make_shared<SimpleTrace>(ch, 0, z));
        gauss_traces.push_back(std::make_shared<SimpleTrace>(ch, 0, z));
    }
    ITrace::vector all = adc_traces;
    for (auto& t : gauss_traces) all.push_back(t);
    for (int i = 0; i < nchan; ++i) raw_idx.push_back(i);
    for (int i = 0; i < nchan; ++i) gauss_idx.push_back(nchan + i);

    auto sf = std::make_shared<SimpleFrame>(42, 0.0, all, 0.5);
    sf->tag_traces("raw", raw_idx);
    sf->tag_traces("gauss", gauss_idx);
    return sf;
}

// Configure L1SPFilterPD with the production kernel, gauss_filter
// disabled (we don't need smearing for a contract pin), and no anode
// (m_anode==null bypasses plane resolution).  Returns the configurable
// pointer so the caller can call configure() once.
Configuration kernel_cfg(IConfigurable::pointer icfg)
{
    auto cfg = icfg->default_configuration();
    cfg["kernels_file"] = KERNEL;
    cfg["gauss_filter"] = "";
    cfg["anode"] = "";
    return cfg;
}

}  // namespace

TEST_CASE("L1SPFilterPD loads production kernel and processes minimal frame") {
    if (!kernel_available()) {
        MESSAGE("SKIP: " << KERNEL << " not found via WIRECELL_PATH");
        return;
    }
    load_plugins();
    auto icfg = Factory::lookup_tn<IConfigurable>("L1SPFilterPD");
    auto cfg = kernel_cfg(icfg);
    icfg->configure(cfg);

    auto ff = Factory::lookup_tn<IFrameFilter>("L1SPFilterPD");
    auto in = make_minimal_frame();
    IFrame::pointer out;
    REQUIRE_NOTHROW(ff->operator()(in, out));
    REQUIRE(out);

    // Output must carry the L1SP output tag (default "l1sp") even when
    // no ROIs trigger, because the operator emits the unmodified traces
    // under that tag.
    auto trace_tags = out->trace_tags();
    bool has_outtag = false;
    for (auto const& t : trace_tags) {
        if (t == "l1sp") { has_outtag = true; break; }
    }
    CHECK(has_outtag);
}

TEST_CASE("L1SPFilterPD dump-mode emits documented NPZ schema") {
    if (!kernel_available()) {
        MESSAGE("SKIP: " << KERNEL << " not found via WIRECELL_PATH");
        return;
    }
    load_plugins();

    // /home/xqian/tmp is the user-blessed scratch area (memory:
    // feedback_tmp_directory).
    const std::string tmp = "/home/xqian/tmp/wct_l1sp_dump_test";
    std::error_code ec;
    std::filesystem::remove_all(tmp, ec);

    auto icfg = Factory::lookup_tn<IConfigurable>("L1SPFilterPD");
    auto cfg = kernel_cfg(icfg);
    cfg["dump_mode"] = true;
    cfg["dump_path"] = tmp;
    cfg["dump_tag"]  = "test";
    icfg->configure(cfg);

    auto ff = Factory::lookup_tn<IFrameFilter>("L1SPFilterPD");
    auto in = make_minimal_frame();
    IFrame::pointer out;
    REQUIRE_NOTHROW(ff->operator()(in, out));

    // Locate the produced NPZ.  Filename pattern is
    //   <tag>_<call_count:04d>_<frame_ident>.npz
    // m_count accumulates across TEST_CASEs because L1SPFilterPD is a
    // factory singleton, so we glob rather than hard-coding 0000.
    std::string found;
    for (auto const& e : std::filesystem::directory_iterator(tmp)) {
        if (e.path().extension() == ".npz") { found = e.path().string(); break; }
    }
    REQUIRE_FALSE(found.empty());

    auto npz = cnpy::npz_load(found);
    auto has = [&](const std::string& key) { return npz.count(key) > 0; };

    // Per-frame scalars (always written)
    CHECK(has("frame_ident"));
    CHECK(has("frame_time"));
    CHECK(has("call_count"));
    CHECK(has("n_rois"));

    // Per-ROI core columns (vectors; may be size 0 on no-trigger input)
    CHECK(has("channel"));
    CHECK(has("roi_start"));
    CHECK(has("roi_end"));
    CHECK(has("nbin_fit"));
    CHECK(has("temp_sum"));
    CHECK(has("temp1_sum"));
    CHECK(has("temp2_sum"));
    CHECK(has("max_val"));
    CHECK(has("min_val"));
    CHECK(has("prev_roi_end"));
    CHECK(has("next_roi_start"));
    CHECK(has("prev_gap"));
    CHECK(has("next_gap"));

    // Tier 1
    CHECK(has("flag"));
    CHECK(has("ratio"));
    CHECK(has("temp_sum_pos"));
    CHECK(has("temp_sum_neg"));
    CHECK(has("n_above_pos"));
    CHECK(has("n_above_neg"));
    // Tier 2
    CHECK(has("argmax_tick"));
    CHECK(has("argmin_tick"));
    CHECK(has("sig_peak"));
    CHECK(has("sig_integral"));
    // Tier 3 — features driving the trigger gate
    CHECK(has("gmax"));
    CHECK(has("gauss_fill"));
    CHECK(has("gauss_fwhm_frac"));
    CHECK(has("roi_energy_frac"));
    CHECK(has("raw_asym_wide"));
    CHECK(has("core_lo"));
    CHECK(has("core_hi"));
    CHECK(has("core_length"));
    CHECK(has("core_fill"));
    CHECK(has("core_fwhm_frac"));
    CHECK(has("core_raw_asym_wide"));
    CHECK(has("flag_l1"));
    // Adjacency
    CHECK(has("flag_l1_adj"));
    CHECK(has("adj_donor_ch"));

    // n_rois on this minimal flat input is 0
    auto n_rois = npz.at("n_rois").as_vec<int32_t>();
    REQUIRE(n_rois.size() == 1);
    CHECK(n_rois[0] == 0);

    std::filesystem::remove_all(tmp, ec);
}
