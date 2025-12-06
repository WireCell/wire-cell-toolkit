
#include "WireCellImg/LCBlobRemoval.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/Projection2D.h"

#include "WireCellAux/SimpleCluster.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Rectangles.h"

#include <iterator>
#include <chrono>

WIRECELL_FACTORY(LCBlobRemoval, WireCell::Img::LCBlobRemoval,
                 WireCell::INamed,
                 WireCell::IClusterFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Img::CS;
using namespace WireCell::Img::Projection2D;

Img::LCBlobRemoval::LCBlobRemoval()
    : Aux::Logger("LCBlobRemoval", "img")
{
}

Img::LCBlobRemoval::~LCBlobRemoval() {}

WireCell::Configuration Img::LCBlobRemoval::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["blob_value_threshold"] = m_blob_thresh.value();
    cfg["blob_error_threshold"] = m_blob_thresh.uncertainty();

    return cfg;
}

void Img::LCBlobRemoval::configure(const WireCell::Configuration& cfg)
{
    float bt_val = get(cfg,"blob_value_threshold", m_blob_thresh.value());
    float bt_err = get(cfg,"blob_error_threshold", m_blob_thresh.uncertainty());
    m_blob_thresh = value_t(bt_val, bt_err);
}

static void dump_cg(const cluster_graph_t& cg, Log::logptr_t& log)
{
    size_t mcount{0}, bcount{0};
    value_t bval;
    for (const auto& vtx : mir(boost::vertices(cg))) {
        const auto& node = cg[vtx];
        if (node.code() == 'b') {
            const auto iblob = get<blob_t>(node.ptr);
            bval += value_t(iblob->value(), iblob->uncertainty());
            ++bcount;
            continue;
        }
        if (node.code() == 'm') {
            const auto imeas = get<meas_t>(node.ptr);
            ++mcount;
        }
    }
    log->debug("cluster graph: vertices={} edges={} #blob={} bval={} #meas={}",
               boost::num_vertices(cg), boost::num_edges(cg),
               bcount, bval, mcount);
}

namespace {
    cluster_graph_t prune(const cluster_graph_t& cg, float threshold)
    {
        cluster_graph_t cg_out;
        
        size_t nblobs = 0;
        std::unordered_map<cluster_vertex_t, cluster_vertex_t> old2new;
        for (const auto& vtx : mir(boost::vertices(cg))) {
            const auto& node = cg[vtx];
            if (node.code() == 'b') {
                const auto iblob = get<blob_t>(node.ptr);
                auto bval = iblob->value();
                if (bval < threshold)
                continue;
            }
            ++nblobs;
            old2new[vtx] = boost::add_vertex(node, cg_out);
        }
        
        if (!nblobs) {
            return cg_out;
        }

        for (auto edge : mir(boost::edges(cg))) {
            auto old_tail = boost::source(edge, cg);
            auto old_head = boost::target(edge, cg);

            auto old_tit = old2new.find(old_tail);
            if (old_tit == old2new.end()) {
                continue;
            }
            auto old_hit = old2new.find(old_head);
            if (old_hit == old2new.end()) {
                continue;
            }
            boost::add_edge(old_tit->second, old_hit->second, cg_out);
        }
        return cg_out;
    }
}

bool Img::LCBlobRemoval::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }
    const auto in_graph = in->graph();
    dump_cg(in_graph, log);
    auto out_graph = prune(in_graph, m_blob_thresh.value());
    dump_cg(out_graph, log);
    out = std::make_shared<Aux::SimpleCluster>(out_graph, in->ident());
    return true;
}


