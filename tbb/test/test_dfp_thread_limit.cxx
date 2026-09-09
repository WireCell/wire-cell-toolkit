// Regression test for the TbbFlow node-parallelism thread limit.
//
// Background (spng-fra.12): DataFlowGraph originally limited node parallelism
// with tbb::global_control(max_allowed_parallelism = max_threads).  That is a
// PROCESS-GLOBAL throttle: it keeps a full hardware-sized TBB worker pool that
// busy-spins while a long node runs (burning every core), AND it composes by
// minimum with any host framework's global_control -- so an embedded WCT would
// clamp the whole host (e.g. Phlex) to WCT's limit.
//
// The fix runs the graph inside a sized tbb::task_arena(max_threads) instead.
// A sized arena is LOCAL and composable: it bounds only WCT's own graph and
// creates no idle workers to spin.  When no limit is set (max_threads=0), the
// graph runs in whatever arena is current on the calling thread -- i.e. it
// INHERITS the host framework's arena.
//
// This test asserts both properties deterministically by observing the arena
// concurrency seen from inside a running node:
//   1. max_threads=N  -> node runs in an arena of concurrency N.
//   2. max_threads=0 inside a host task_arena(M) -> node runs at concurrency M
//      (WCT defers to the host).
//
// Under the old global_control implementation, case 1 would observe the default
// (hardware) concurrency, not N, so this test would fail -- as intended.

#include "tbb_mock.h"

#include "WireCellTbb/DataFlowGraph.h"
#include "WireCellUtil/Testing.h"

#include <tbb/task_arena.h>
#include <tbb/info.h>

#include <atomic>
#include <iostream>

using namespace WireCell;

// A depo sink that records the TBB arena concurrency it executes under.
static std::atomic<int> g_observed{-1};

class ConcurrencySink : public IDepoSink, public Aux::Logger {
  public:
    ConcurrencySink() : Aux::Logger("ConcurrencySink", "test") {}
    virtual ~ConcurrencySink() {}
    virtual bool operator()(const input_pointer& /*depo*/)
    {
        g_observed.store(tbb::this_task_arena::max_concurrency());
        return true;
    }
};

// Build source -> sink, run with the given max_threads, return observed arena
// concurrency.  A fresh graph/nodes per call (connect and the source are
// one-shot).
static int run_and_observe(int max_threads)
{
    g_observed.store(-1);
    WireCellTbb::DataFlowGraph dfg;
    Configuration cfg = dfg.default_configuration();
    cfg["max_threads"] = max_threads;
    cfg["summary"] = 0;                 // keep the test quiet
    dfg.configure(cfg);

    INode::pointer source(new WireCellTbb::MockDepoSource(3));
    INode::pointer sink(new ConcurrencySink);
    Assert(dfg.connect(source, sink));
    Assert(dfg.run());
    return g_observed.load();
}

int main()
{
    const int hw = tbb::info::default_concurrency();
    std::cerr << "default_concurrency=" << hw << std::endl;

    // Case 1: an explicit limit runs the graph in a sized arena.
    const int limit = 2;
    int obs1 = run_and_observe(limit);
    std::cerr << "max_threads=" << limit << " -> observed concurrency " << obs1 << std::endl;
    Assert(obs1 == limit);

    // Case 2 (embedding): with no WCT limit, run inside a host arena and confirm
    // WCT inherits the host's concurrency rather than overriding it.  Only
    // meaningful when the host size differs from the machine default.
    const int host = (hw >= 3) ? 3 : 1;
    int obs2 = -1;
    tbb::task_arena host_arena(host);
    host_arena.execute([&]() { obs2 = run_and_observe(0); });
    std::cerr << "host arena=" << host << ", max_threads=0 -> observed concurrency "
              << obs2 << std::endl;
    Assert(obs2 == host);

    std::cerr << "test_dfp_thread_limit OK\n";
    return 0;
}
