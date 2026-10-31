// Test that Graph::execute() (used by Pgrapher) is a faithful, deadlock-free
// reimplementation of the original Pgrapher scheduler, kept here as
// reference_execute(): for the same graph it must fire the identical sequence
// of nodes and deliver the identical data to every sink -- only the scheduling
// algorithm changed (O((V+E) log V) vs O(V^2)), not the execution policy.

#include "WireCellPgraph/Graph.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <boost/any.hpp>
#include <memory>
#include <vector>

using namespace WireCell;
using namespace WireCell::Pgraph;

namespace {

    // Shared log of node ids, appended to only on a SUCCESSFUL firing.
    struct Rec {
        std::vector<int> fires;
    };

    struct TestNode : public Node {
        int m_id;
        Rec* m_rec;
        TestNode(int id, Rec* rec) : m_id(id), m_rec(rec) {}
        std::string ident() const override { return instance_name(); }
        std::string instance_name() const override { return "n" + std::to_string(m_id); }
        std::string cpptype_name() const override { return "TestNode"; }
        void fired() { m_rec->fires.push_back(m_id); }
    };

    // Emits 0..n-1 then goes dry.  Self-guards against overproduction so no
    // more than one datum ever waits on its output edge (the in-flight policy).
    struct Src : public TestNode {
        int m_num{0}, m_end;
        Src(int id, int n, Rec* rec) : TestNode(id, rec), m_end(n)
        {
            m_ports[Port::output].push_back(Port(this, Port::output, "int"));
        }
        bool operator()() override
        {
            if (oport().size()) return false;   // output still pending
            if (m_num >= m_end) return false;   // exhausted
            Data d = m_num++;
            oport().put(d);
            fired();
            return true;
        }
    };

    struct Func : public TestNode {
        Func(int id, Rec* rec) : TestNode(id, rec)
        {
            m_ports[Port::input].push_back(Port(this, Port::input, "int"));
            m_ports[Port::output].push_back(Port(this, Port::output, "int"));
        }
        bool operator()() override
        {
            if (oport().size()) return false;
            if (iport().empty()) return false;
            Data d = iport().get();
            oport().put(d);
            fired();
            return true;
        }
    };

    struct Fanout : public TestNode {
        Fanout(int id, size_t nout, Rec* rec) : TestNode(id, rec)
        {
            m_ports[Port::input].push_back(Port(this, Port::input, "int"));
            for (size_t ind = 0; ind < nout; ++ind) {
                m_ports[Port::output].push_back(Port(this, Port::output, "int"));
            }
        }
        bool operator()() override
        {
            if (iport().empty()) return false;
            bool full = true;
            for (auto& p : output_ports()) {
                if (p.empty()) { full = false; break; }
            }
            if (full) return false;             // all outputs still pending
            auto in = iport().get();
            for (auto& p : output_ports()) {
                Data d = in;
                p.put(d);
            }
            fired();
            return true;
        }
    };

    struct Sink : public TestNode {
        std::vector<int>* m_got;
        Sink(int id, std::vector<int>* got, Rec* rec) : TestNode(id, rec), m_got(got)
        {
            m_ports[Port::input].push_back(Port(this, Port::input, "int"));
        }
        bool operator()() override
        {
            if (iport().empty()) return false;
            m_got->push_back(boost::any_cast<int>(iport().get()));
            fired();
            return true;
        }
    };

    // The original Graph::execute(): reverse-scan the topological order,
    // fire the first (sink-most) runnable node, restart from the sinks.
    void reference_execute(Graph& g)
    {
        auto nodes = g.sort_kahn();
        bool did_something = true;
        while (did_something) {
            did_something = false;
            for (auto nit = nodes.rbegin(); nit != nodes.rend(); ++nit) {
                if (g.call_node(*nit)) {
                    did_something = true;
                    break;
                }
            }
        }
    }

    // Build source -> func -> fanout(W) -> W sinks, run with either scheduler,
    // and return the firing log plus each sink's received values.
    struct Result {
        std::vector<int> fires;
        std::vector<std::vector<int>> got;
    };

    Result run(size_t W, int nvals, bool reference)
    {
        Rec rec;
        std::vector<std::vector<int>> got(W);

        int id = 0;
        auto src = std::make_shared<Src>(id++, nvals, &rec);
        auto fun = std::make_shared<Func>(id++, &rec);
        auto fan = std::make_shared<Fanout>(id++, W, &rec);
        std::vector<std::shared_ptr<Sink>> sinks;
        for (size_t ind = 0; ind < W; ++ind) {
            sinks.push_back(std::make_shared<Sink>(id++, &got[ind], &rec));
        }

        Graph g;
        g.connect(src.get(), fun.get());
        g.connect(fun.get(), fan.get());
        for (size_t ind = 0; ind < W; ++ind) {
            g.connect(fan.get(), sinks[ind].get(), ind, 0);
        }

        if (reference) {
            reference_execute(g);
        }
        else {
            g.execute();
        }

        return Result{rec.fires, got};
    }

}  // namespace

TEST_CASE("pgraph execute matches reference on a fanout tree")
{
    const size_t W = 7;
    const int nvals = 5;

    auto a = run(W, nvals, /*reference=*/true);   // original rescan scheduler
    auto b = run(W, nvals, /*reference=*/false);  // Graph::execute()

    // Identical sequence of successful node firings: same policy, same order.
    CHECK(a.fires == b.fires);

    // Every sink received every value, in order, under both engines.
    REQUIRE(a.got.size() == W);
    REQUIRE(b.got.size() == W);
    std::vector<int> expect;
    for (int v = 0; v < nvals; ++v) expect.push_back(v);
    for (size_t ind = 0; ind < W; ++ind) {
        CHECK(a.got[ind] == expect);
        CHECK(b.got[ind] == expect);
    }
}

TEST_CASE("pgraph execute terminates and matches reference for several fan widths")
{
    for (size_t W : {1u, 2u, 3u, 10u, 25u}) {
        auto a = run(W, 3, true);
        auto b = run(W, 3, false);
        CHECK(a.fires == b.fires);
        CHECK(a.got == b.got);
    }
}
