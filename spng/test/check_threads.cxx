// check_threads - explore torch/MKL thread control and core-time measurement,
// especially when a torch op runs on a TBB worker thread (as under TbbFlow).
//
// Background: under the TbbFlow engine, SPNG wire-cell jobs show a per-node
// core-sec/wall-sec ratio of ~= (#hardware threads) for EVERY node, including
// non-parallel ones (file sinks), while OSP jobs on the same engine do not.
// The shell's own `time` confirms the process really burns ~64x CPU.  This toy
// isolates the phenomenon so we can (a) find a control that actually caps the
// CPU torch uses on a worker thread and (b) find a core-time metric that is not
// fooled by idle-thread spinning.
//
// Build: this is a check_*.cxx program -> built (not run) as WireCellSpng__check_threads.
// Run examples:
//   ./WireCellSpng__check_threads --op matmul --where main
//   ./WireCellSpng__check_threads --op matmul --where tbb --tbb-max 1
//   ./WireCellSpng__check_threads --op matmul --where tbb --tbb-max 1 --omp 1 --apply-in-op
//   ./WireCellSpng__check_threads --op matmul --where tbb --tbb-max 1 --idle-ms 500
//
// The metrics printed per phase:
//   wall   : steady_clock elapsed seconds
//   clock  : std::clock() delta / CLOCKS_PER_SEC  (process CPU; what pgraph/tbb log)
//   ru_self: getrusage(RUSAGE_SELF) user+sys delta (process CPU, another route)
//   ru_thr : getrusage(RUSAGE_THREAD) user+sys delta on the CALLING thread only
// and the thread-count views at/inside the op: at::get_num_threads(),
// omp_get_max_threads(), mkl_get_max_threads().

#include "WireCellSpng/Torch.h"        // sole torch include site for SPNG
#include "WireCellSpng/TorchScript.h"  // torch::jit::load / <torch/script.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_TBB
#include <tbb/global_control.h>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/flow_graph.h>
#endif

#include <sys/resource.h>
#include <dlfcn.h>

#include <chrono>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <string>
#include <thread>
#include <iostream>
#include <functional>

// ------- MKL controls via dlsym (avoid an MKL header/link dependency) -------
namespace mkl {
    using set_t = void (*)(int);
    using setlocal_t = int (*)(int);
    using get_t = int (*)(void);

    static set_t      set_num_threads       = (set_t)      dlsym(RTLD_DEFAULT, "mkl_set_num_threads");
    static setlocal_t set_num_threads_local = (setlocal_t) dlsym(RTLD_DEFAULT, "mkl_set_num_threads_local");
    static get_t      get_max_threads       = (get_t)      dlsym(RTLD_DEFAULT, "mkl_get_max_threads");

    int max() { return get_max_threads ? get_max_threads() : -1; }
}

// ------- timing -------
struct Sample {
    std::chrono::steady_clock::time_point wall;
    std::clock_t clk;
    double ru_self;   // seconds
    double ru_thr;    // seconds

    static double rusage_secs(int who) {
        struct rusage ru;
        getrusage(who, &ru);
        double u = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6;
        double s = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
        return u + s;
    }
    static Sample now() {
        Sample s;
        s.wall = std::chrono::steady_clock::now();
        s.clk = std::clock();
        s.ru_self = rusage_secs(RUSAGE_SELF);
#ifdef RUSAGE_THREAD
        s.ru_thr = rusage_secs(RUSAGE_THREAD);
#else
        s.ru_thr = 0;
#endif
        return s;
    }
};

static void report(const std::string& tag, const Sample& a, const Sample& b)
{
    double wall = std::chrono::duration<double>(b.wall - a.wall).count();
    double clk  = double(b.clk - a.clk) / CLOCKS_PER_SEC;
    double self = b.ru_self - a.ru_self;
    double thr  = b.ru_thr  - a.ru_thr;
    auto r = [&](double c){ return wall > 1e-9 ? c / wall : 0.0; };
    std::printf("  %-10s wall=%8.4f  clock=%9.4f (%.1fx)  ru_self=%9.4f (%.1fx)  ru_thr=%8.4f (%.1fx)\n",
                tag.c_str(), wall, clk, r(clk), self, r(self), thr, r(thr));
}

// ------- thread control -------
struct Controls {
    int torch = 0;   // at::set_num_threads (intra-op)
    int interop = 0; // at::set_num_interop_threads (once, before any op)
    int omp = 0;     // omp_set_num_threads
    int mkl = 0;     // mkl_set_num_threads (global)
    int mkl_local = 0; // mkl_set_num_threads_local (per-thread)

    void apply(const char* where) const {
        if (torch > 0) at::set_num_threads(torch);
#ifdef _OPENMP
        if (omp > 0) omp_set_num_threads(omp);
#endif
        if (mkl > 0 && mkl::set_num_threads) mkl::set_num_threads(mkl);
        if (mkl_local > 0 && mkl::set_num_threads_local) mkl::set_num_threads_local(mkl_local);
        (void) where;
    }
};

static void print_thread_view(const char* tag)
{
    int omp_max = -1;
#ifdef _OPENMP
    omp_max = omp_get_max_threads();
#endif
    std::printf("  [%s] at::get_num_threads=%d  at::get_num_interop=%d  omp_get_max_threads=%d  mkl_get_max_threads=%d\n",
                tag, (int)at::get_num_threads(), (int)at::get_num_interop_threads(), omp_max, mkl::max());
}

// ------- the compute op -------
struct Op {
    std::string kind = "matmul";  // matmul | conv | jit  (jit = scripted conv module)
    int64_t N = 2048;     // matmul square dim
    // conv params (roiuniter-ish: 1 channel, big image)
    int64_t C = 8, H = 1500, W = 800, K = 3;

    torch::Tensor A, B, img, wgt;
    mutable torch::jit::Module mod;
    double busy_s = 2.0;   // seconds for the "busy" (guaranteed single-thread) op
    std::string model_path;                 // for kind=="model"
    std::vector<int64_t> in_shape{1,1,800,1500};

    void init() {
        torch::manual_seed(1);
        if (kind == "busy") return;
        if (kind == "model") {
            mod = torch::jit::load(model_path);
            mod.eval();
            img = torch::randn(in_shape);
            return;
        }
        if (kind == "matmul") {
            A = torch::randn({N, N});
            B = torch::randn({N, N});
        } else {
            img = torch::randn({1, C, H, W});
            wgt = torch::randn({C, C, K, K});
        }
        if (kind == "jit") {
            // Build a TorchScript module doing the same conv, to mirror how the
            // SPNG roiuniter forward runs (torch::jit::Module::forward), which
            // may use different parallelism than eager ops.
            mod = torch::jit::Module("m");
            mod.register_parameter("w", wgt, /*is_buffer*/false);
            mod.define(R"JIT(
def forward(self, x):
    return torch.conv2d(x, self.w, None, [1, 1], [1, 1])
)JIT");
        }
    }
    // Run once; returns something to prevent optimizing away.
    double run() const {
        if (kind == "busy") {
            // Guaranteed single-threaded CPU work for busy_s seconds.  Isolates
            // the execution engine's threading from torch's own.
            auto t0 = std::chrono::steady_clock::now();
            volatile double x = 0;
            while (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() < busy_s) {
                for (int i = 0; i < 200000; ++i) x += i * 1e-9;
            }
            return x;
        }
        torch::NoGradGuard ng;
        torch::Tensor out;
        if (kind == "matmul") out = torch::matmul(A, B);
        else if (kind == "conv") out = torch::conv2d(img, wgt, /*bias*/{}, /*stride*/1, /*padding*/1);
        else /* jit or model */ out = mod.forward({img}).toTensor();
        return out.sum().item<double>();
    }
};

// ------- CLI -------
static int argi(int argc, char** argv, int& i, const char* name) {
    if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", name); std::exit(2); }
    return std::atoi(argv[++i]);
}

int main(int argc, char** argv)
{
    Op op;
    Controls ctl;
    std::string where = "main";     // main | thread | tbb
    int tbb_max = 0;                // global_control max_allowed_parallelism
    int arena = 0;                  // task_arena concurrency (0 = default)
    bool apply_in_op = false;       // apply controls on the executing thread inside the op
    int iters = 3, warmup = 1;
    int idle_ms = 0;
    int fgpar = 0;                  // >0: flowgraph with this many parallel nodes
    int items = 0;                  // >0: fgsrc mode, input_node produces this many items
    std::string limit_mode = "gc";  // gc (global_control) | arena (task_arena)

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a == "--op")        op.kind = argv[++i];
        else if (a == "--busy-s")    op.busy_s = std::atof(argv[++i]);
        else if (a == "--model")     op.model_path = argv[++i];
        else if (a == "--in-shape") {
            op.in_shape.clear();
            std::string s = argv[++i], tok;
            for (char c : s) {
                if (c == ',') { op.in_shape.push_back(std::atoll(tok.c_str())); tok.clear(); }
                else tok += c;
            }
            if (!tok.empty()) op.in_shape.push_back(std::atoll(tok.c_str()));
        }
        else if (a == "--N")         op.N = argi(argc, argv, i, "--N");
        else if (a == "--iters")     iters = argi(argc, argv, i, "--iters");
        else if (a == "--warmup")    warmup = argi(argc, argv, i, "--warmup");
        else if (a == "--torch")     ctl.torch = argi(argc, argv, i, "--torch");
        else if (a == "--interop")   ctl.interop = argi(argc, argv, i, "--interop");
        else if (a == "--omp")       ctl.omp = argi(argc, argv, i, "--omp");
        else if (a == "--mkl")       ctl.mkl = argi(argc, argv, i, "--mkl");
        else if (a == "--mkl-local") ctl.mkl_local = argi(argc, argv, i, "--mkl-local");
        else if (a == "--where")     where = argv[++i];
        else if (a == "--tbb-max")   tbb_max = argi(argc, argv, i, "--tbb-max");
        else if (a == "--arena")     arena = argi(argc, argv, i, "--arena");
        else if (a == "--apply-in-op") apply_in_op = true;
        else if (a == "--idle-ms")   idle_ms = argi(argc, argv, i, "--idle-ms");
        else if (a == "--fgpar")     fgpar = argi(argc, argv, i, "--fgpar");
        else if (a == "--items")     items = argi(argc, argv, i, "--items");
        else if (a == "--limit-mode") limit_mode = argv[++i];
        else { std::fprintf(stderr, "unknown arg: %s\n", a.c_str()); return 2; }
    }

    std::printf("config: op=%s where=%s tbb_max=%d arena=%d apply_in_op=%d "
                "torch=%d omp=%d mkl=%d mkl_local=%d iters=%d warmup=%d idle_ms=%d\n",
                op.kind.c_str(), where.c_str(), tbb_max, arena, (int)apply_in_op,
                ctl.torch, ctl.omp, ctl.mkl, ctl.mkl_local, iters, warmup, idle_ms);

    // Inter-op threads can only be set ONCE, before any inter-op work.  Do it
    // first, regardless of apply_in_op.
    if (ctl.interop > 0) {
        try { at::set_num_interop_threads(ctl.interop); }
        catch (const std::exception& e) {
            std::fprintf(stderr, "set_num_interop_threads failed: %s\n", e.what());
        }
    }

    // Apply the rest up front (on the main thread) unless deferred to the op.
    if (!apply_in_op) ctl.apply("main-preop");
    print_thread_view("startup");

    op.init();

#ifdef HAVE_TBB
    // Concurrency limiting: "gc" mimics TbbFlow today (global_control on
    // max_allowed_parallelism -- excluded workers may busy-spin when tasks are
    // pending); "arena" limits the WORKER COUNT via a sized task_arena so no
    // excess workers exist to spin.
    std::unique_ptr<tbb::global_control> gc;
    if (tbb_max > 0 && limit_mode == "gc") {
        gc.reset(new tbb::global_control(tbb::global_control::max_allowed_parallelism, tbb_max));
    }
#endif

    double sink = 0;

    // The op body, possibly re-applying controls on the executing thread.
    auto body = [&]() {
        if (apply_in_op) {
            ctl.apply("in-op");
            print_thread_view("in-op-thread");
        }
        sink += op.run();
    };

    // A runner that places `body` on the requested execution context.
    auto run_once = [&]() {
        if (where == "main") {
            body();
        }
        else if (where == "thread") {
            std::thread t(body);
            t.join();
        }
        else if (where == "tbb") {
#ifdef HAVE_TBB
            tbb::task_arena ta(arena > 0 ? arena : tbb::task_arena::automatic);
            ta.execute([&]() {
                // parallel_for with one item forces execution on a TBB worker.
                tbb::parallel_for(0, 1, [&](int) { body(); });
            });
#else
            std::fprintf(stderr, "built without TBB; --where tbb unavailable\n");
            std::exit(3);
#endif
        }
        else if (where == "flowgraph") {
#ifdef HAVE_TBB
            // Mimic TbbFlow: a tbb::flow::graph with one function_node running the
            // op, driven to completion by wait_for_all().  While the graph is
            // active, do the other TBB workers spin-wait and burn CPU?
            tbb::flow::graph g;
            tbb::flow::function_node<int, int> node(
                g, tbb::flow::serial, [&](int) -> int { body(); return 0; });
            node.try_put(0);
            g.wait_for_all();
#else
            std::fprintf(stderr, "built without TBB; --where flowgraph unavailable\n");
            std::exit(3);
#endif
        }
        else { std::fprintf(stderr, "bad --where %s\n", where.c_str()); std::exit(2); }
    };

    // Multi-node flow graph: a broadcast feeding `fgpar` parallel function
    // nodes; node 0 runs the long op, the rest run a short op.  This creates
    // several ready tasks at once so a concurrency limit of 1 leaves the other
    // workers with pending-but-unrunnable work -- the condition under which TBB
    // workers busy-spin.  Optionally run inside a sized task_arena.
    auto run_fgpar = [&]() {
#ifdef HAVE_TBB
        auto build_and_run = [&]() {
            tbb::flow::graph g;
            tbb::flow::broadcast_node<int> bcast(g);
            std::vector<std::unique_ptr<tbb::flow::function_node<int,int>>> nodes;
            for (int n = 0; n < fgpar; ++n) {
                // unlimited concurrency: all fgpar nodes become runnable at once
                // on broadcast, so a throttle of 1 leaves fgpar-1 tasks pending
                // but runnable -- the condition that makes TBB workers busy-spin.
                nodes.emplace_back(new tbb::flow::function_node<int,int>(
                    g, tbb::flow::unlimited, [&](int) -> int { sink += op.run(); return 0; }));
                tbb::flow::make_edge(bcast, *nodes.back());
            }
            bcast.try_put(0);
            g.wait_for_all();
        };
        if (tbb_max > 0 && limit_mode == "arena") {
            tbb::task_arena ta(tbb_max);
            ta.execute(build_and_run);
        } else {
            build_and_run();
        }
#endif
    };

    // fgsrc: mimic WCT's TbbFlow more faithfully -- a tbb::flow::input_node
    // source (as NodeWrapper.h uses) producing `items` messages into a slow
    // serial function_node, driven by wait_for_all, limited by gc or arena.
    auto run_fgsrc = [&]() {
#ifdef HAVE_TBB
        auto build_and_run = [&]() {
            tbb::flow::graph g;
            int produced = 0;
            tbb::flow::input_node<int> src(g, [&](tbb::flow_control& fc) -> int {
                if (produced >= items) { fc.stop(); return 0; }
                return produced++;
            });
            tbb::flow::function_node<int,int> slow(
                g, tbb::flow::serial, [&](int) -> int {
                    std::printf("  [fgsrc-node] this_task_arena::max_concurrency=%d\n",
                                tbb::this_task_arena::max_concurrency());
                    sink += op.run(); return 0; });
            tbb::flow::make_edge(src, slow);
            src.activate();
            g.wait_for_all();
        };
        if (tbb_max > 0 && limit_mode == "arena") {
            tbb::task_arena ta(tbb_max);
            ta.execute(build_and_run);
        } else {
            build_and_run();
        }
#endif
    };

    auto do_one = [&]() {
        if (items > 0) run_fgsrc();
        else if (fgpar > 0) run_fgpar();
        else run_once();
    };

    for (int w = 0; w < warmup; ++w) do_one();

    std::printf("timed (%d iters):\n", iters);
    Sample a = Sample::now();
    for (int k = 0; k < iters; ++k) do_one();
    Sample b = Sample::now();
    report("op", a, b);

    // Idle probe: with the op's thread pool now warm, measure CPU burned while
    // the process does nothing.  Nonzero ru_self/clock here == background spin.
    if (idle_ms > 0) {
        Sample c = Sample::now();
        std::this_thread::sleep_for(std::chrono::milliseconds(idle_ms));
        Sample d = Sample::now();
        report("idle", c, d);
    }

    print_thread_view("end");
    std::printf("sink=%g\n", sink);
    return 0;
}
