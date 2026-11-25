// Pgrapher runs a hydra node when ANY of its input ports holds data.
//
// A hydra (IHydraNode) may consume and buffer its input (eg Aux::FrameSync,
// which must, else the leftover input is handed back and duplicated).  If the
// wrapper waited for ALL input ports to be non-empty, a datum arriving on one
// port after the hydra emptied the others would never be delivered: here the
// frame and its EOS come on port 0 while port 1 carries only EOS, and the
// merged EOS never reached the sink.  In a job that runs the graph once per
// event (eg wire-cell-phlex) the stranded EOS then fires on the next run
// before that run's frame arrives, and the frame is lost.  TbbFlow already
// calls a hydra per message on any port.

#include "WireCellPgraph/Graph.h"
#include "WireCellPgraph/Wrappers.h"
#include "WireCellAux/FrameSync.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/IFrameSource.h"
#include "WireCellUtil/doctest.h"

#include <deque>
#include <memory>
#include <vector>

using namespace WireCell;

namespace {

    // Emits the given frames, then EOS, then nothing.
    struct FrameList : public IFrameSource {
        std::deque<IFrame::pointer> m_frames;
        bool m_eos_sent{false};
        explicit FrameList(std::vector<int> idents)
        {
            for (int ident : idents) {
                m_frames.push_back(std::make_shared<Aux::SimpleFrame>(ident));
            }
        }
        bool operator()(IFrame::pointer& out) override
        {
            out = nullptr;
            if (!m_frames.empty()) {
                out = m_frames.front();
                m_frames.pop_front();
                return true;
            }
            if (!m_eos_sent) {
                m_eos_sent = true;
                return true;
            }
            return false;
        }
    };

    // Records frame idents, -1 for EOS.
    struct Collect : public IFrameSink {
        std::vector<int> got;
        bool operator()(const IFrame::pointer& in) override
        {
            got.push_back(in ? in->ident() : -1);
            return true;
        }
    };

}  // namespace

TEST_CASE("pgraph hydra runs when any input port has data")
{
    auto port0 = std::make_shared<FrameList>(std::vector<int>{1});
    auto port1 = std::make_shared<FrameList>(std::vector<int>{});
    auto sync = std::make_shared<Aux::FrameSync>(2);
    auto sink = std::make_shared<Collect>();

    Pgraph::Source n0(port0), n1(port1);
    Pgraph::Hydra hy(sync);
    Pgraph::Sink ns(sink);

    Pgraph::Graph graph;
    REQUIRE(graph.connect(&n0, &hy, 0, 0));
    REQUIRE(graph.connect(&n1, &hy, 0, 1));
    REQUIRE(graph.connect(&hy, &ns, 0, 0));
    graph.execute();

    // The frame and then the merged EOS, both delivered in this one run.
    CHECK(sink->got == std::vector<int>{1, -1});
}
