#include "WireCellImg/NaiveStriper.h"
#include "WireCellImg/ImgData.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Graph.h"



WIRECELL_FACTORY(NaiveStriper, WireCell::Img::NaiveStriper, WireCell::ISliceStriper, WireCell::IConfigurable)

using namespace WireCell;

Img::NaiveStriper::~NaiveStriper() {}

WireCell::Configuration Img::NaiveStriper::default_configuration() const
{
    Configuration cfg;

    // The number of intersticial wires which must exist without
    // activity to consider two wires to be non-adjacent.
    cfg["gap"] = 1;

    return cfg;
}

void Img::NaiveStriper::configure(const WireCell::Configuration& cfg) { m_gap = get(cfg, "gap", 1); }

bool Img::NaiveStriper::operator()(const input_pointer& slice, output_pointer& out)
{
    out = nullptr;
    if (!slice) {
        return true;  // eos
    }

    // The graph connects channels to attached wires and wires to
    // their adjacent neighbor in the plane and along the positive
    // pitch direction.

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_t;
    typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
    graph_t graph;

    // Boost graph uses simple numbers as the "node object".  We could
    // attach vertex properties to nodes to hold either the IWire or
    // the IChannel+charge values or we can keep or own lookup tables.
    // The former is probably faster but for now, we these lookups out
    // in the open.
    std::unordered_map<IWire::pointer, IChannel::pointer> wire_to_chan;
    std::unordered_map<IChannel::pointer, vertex_t> chan_to_node;
    std::unordered_map<vertex_t, ISlice::pair_t> node_to_chanval;
    std::unordered_map<int, IWireIndexSet> hit_wires;

    // Fill channel nodes in graph and group wires by plane
    for (const auto& cv : slice->activity()) {
        auto ichan = cv.first;
        auto node = boost::add_vertex(graph);
        chan_to_node[ichan] = node;
        node_to_chanval[node] = cv;
        for (auto iwire : ichan->wires()) {
            const auto pid = iwire->planeid();
            hit_wires[pid.ident()].insert(iwire);
            wire_to_chan[iwire] = ichan;
        }
    }

    // Loop over ordered wires per plane and make edges.
    for (const auto& phw : hit_wires) {
        int last_ind = -1;
        vertex_t last_wire = 0;
        for (auto iwire : phw.second) {
            auto ichan = wire_to_chan[iwire];
            vertex_t chan_node = chan_to_node[ichan];
            vertex_t wire_node = boost::add_vertex(graph);
            const int this_index = iwire->index();

            // no matter what, chan->node.
            boost::add_edge(chan_node, wire_node, graph);

            const size_t dind = this_index - last_ind;
            if (last_ind >= 0 and dind <= m_gap) {
                boost::add_edge(last_wire, wire_node, graph);
            }
            last_ind = this_index;
            last_wire = wire_node;
        }
    }

    // Here's the heavy lifing.  Stripes are understood to be formed as
    // the channels found by looking for the "connected subgraphs".
    // Like the graph itself, this neesds some looksup to translate
    // between Boost Graph's subgraph index and a corresponding stripe.
    std::unordered_map<vertex_t, int> subclusters;
    std::unordered_map<int, Img::Data::Stripe*> cluster_to_stripe;
    boost::connected_components(graph, boost::make_assoc_property_map(subclusters));

    // Collect channels of like cluster number into stripes
    for (auto& p : subclusters) {
        auto ncvit = node_to_chanval.find(p.first);
        if (ncvit == node_to_chanval.end()) {
            continue;
        }
        auto& cv = ncvit->second;
        auto stripe = cluster_to_stripe[p.second];
        if (!stripe) {
            cluster_to_stripe[p.second] = stripe = new Img::Data::Stripe(p.first);
        }
        stripe->append(cv.first, cv.second);
    }

    auto sliceset = new Img::Data::StripeSet(slice->ident());
    for (auto ss : cluster_to_stripe) {
        sliceset->push_back(IStripe::pointer(ss.second));
    }

    out = IStripeSet::pointer(sliceset);
    return true;
}
