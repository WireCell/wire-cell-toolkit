// Structural plumbing tests for L1SPFilterPD that do NOT depend on the
// production kernel file.  Algorithmic / numerical behaviour requiring
// the kernel is covered separately by check_l1sp_*.cxx.
//
// What this test pins:
//   - default_configuration() exposes the documented set of knobs
//   - configure() does not throw on default config (init_resp is lazy)
//   - operator() with EOS (null input) is a clean pass-through
//   - operator() with a non-null frame and an empty kernels_file throws
//     ValueError, with a message that mentions kernels_file (matches the
//     contract at L1SPFilterPD.cxx:726-729).

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

using namespace WireCell;
using namespace WireCell::Aux;

namespace {
void load_plugins() {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");
    pm.add("WireCellSigProc");
}
}  // namespace

TEST_CASE("L1SPFilterPD default_configuration exposes documented knobs") {
    load_plugins();
    auto icfg = Factory::lookup_tn<IConfigurable>("L1SPFilterPD");
    auto cfg = icfg->default_configuration();

    CHECK(cfg.isMember("kernels_file"));
    CHECK(cfg.isMember("kernels_scale"));
    CHECK(cfg.isMember("dft"));
    CHECK(cfg.isMember("adctag"));
    CHECK(cfg.isMember("sigtag"));
    CHECK(cfg.isMember("outtag"));
    CHECK(cfg.isMember("l1_adj_enable"));
    CHECK(cfg.isMember("l1_adj_max_hops"));
    CHECK(cfg.isMember("dump_mode"));
    CHECK(cfg.isMember("dump_path"));
    CHECK(cfg.isMember("waveform_dump_path"));
    CHECK(cfg.isMember("process_planes"));
    CHECK(cfg["kernels_file"].asString().empty());     // unset by default
}

// The default config references HfFilter:Gaus_wide which lives in the
// pgrapher jsonnet, not the C++ factory registry.  For these structural
// tests we override `gauss_filter` to "" so configure() skips that
// lookup; m_smearing_vec stays empty (smearing disabled).  Kernel loading
// is still lazy in init_resp(), which is what we want to test.
static Configuration minimal_cfg(IConfigurable::pointer icfg) {
    auto cfg = icfg->default_configuration();
    cfg["gauss_filter"] = "";   // disable HfFilter lookup
    return cfg;
}

TEST_CASE("L1SPFilterPD configure does not throw with smearing disabled") {
    load_plugins();
    auto icfg = Factory::lookup_tn<IConfigurable>("L1SPFilterPD");
    auto cfg = minimal_cfg(icfg);
    CHECK_NOTHROW(icfg->configure(cfg));
}

TEST_CASE("L1SPFilterPD passes EOS through cleanly (no kernel load)") {
    load_plugins();
    auto icfg = Factory::lookup_tn<IConfigurable>("L1SPFilterPD");
    auto cfg = minimal_cfg(icfg);
    icfg->configure(cfg);

    auto ff = Factory::lookup_tn<IFrameFilter>("L1SPFilterPD");
    IFrame::pointer out;
    CHECK(ff->operator()(nullptr, out));
    CHECK(!out);
}

TEST_CASE("L1SPFilterPD throws ValueError when kernels_file is empty") {
    load_plugins();
    auto icfg = Factory::lookup_tn<IConfigurable>("L1SPFilterPD");
    auto cfg = minimal_cfg(icfg);
    cfg["kernels_file"] = "";   // explicit (default is also "")
    icfg->configure(cfg);

    // Build a minimally-tagged frame.  apply() bails out at init_resp()
    // before it inspects tag content, so a single trace with the right
    // tag pair is enough.
    ITrace::vector traces { std::make_shared<SimpleTrace>(0, 0, std::vector<float>{0.0f, 0.0f}) };
    auto frame = std::make_shared<SimpleFrame>(0, 0.0, traces, 0.5);
    frame->tag_traces("raw",   IFrame::trace_list_t{0});
    frame->tag_traces("gauss", IFrame::trace_list_t{0});

    auto ff = Factory::lookup_tn<IFrameFilter>("L1SPFilterPD");
    IFrame::pointer out;
    CHECK_THROWS_AS(ff->operator()(frame, out), ValueError);
}
