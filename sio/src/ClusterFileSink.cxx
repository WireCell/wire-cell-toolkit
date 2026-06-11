#include "WireCellSio/ClusterFileSink.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/ClusterArrays.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/BlobTools.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/custard/custard_boost.hpp"

WIRECELL_FACTORY(ClusterFileSink, WireCell::Sio::ClusterFileSink,
                 WireCell::INamed,
                 WireCell::IClusterSink, WireCell::ITerminal,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using WireCell::GraphTools::mir;

Sio::ClusterFileSink::ClusterFileSink()
    : Aux::Logger("ClusterFileSink", "io")
{
}

Sio::ClusterFileSink::~ClusterFileSink()
{
}

void Sio::ClusterFileSink::finalize()
{
    log->debug("closing {} after {} calls", m_outname, m_count);
    m_out.pop();
}

WireCell::Configuration Sio::ClusterFileSink::default_configuration() const
{
    Configuration cfg;
    cfg["outname"] = m_outname;
    cfg["prefix"] = m_prefix;
    cfg["format"] = m_format;

    return cfg;
}

std::string Sio::ClusterFileSink::fqn(const ICluster& cluster, std::string name, std::string ext)
{
    std::stringstream ss;
    ss << m_prefix << "_" << cluster.ident() << "_" << name << "." << ext;
    return ss.str();
}

void Sio::ClusterFileSink::jsonify(const ICluster& cluster)
{
    auto top = Aux::jsonify(cluster.graph());
    std::stringstream topss;
    topss << top;
    auto tops = topss.str();
    m_out << "name " << fqn(cluster, "graph", "json") << "\n"
          << "body " << tops.size() << "\n" << tops.data();
    m_out.flush();
}



void Sio::ClusterFileSink::dotify(const ICluster& cluster)
{
    auto top = Aux::dotify(cluster.graph());
    std::stringstream topss;
    topss << top;
    auto tops = topss.str();
    m_out << "name " << fqn(cluster, "graph", "dot") << "\n"
          << "body " << tops.size() << "\n" << tops.data();
    m_out.flush();
}


void Sio::ClusterFileSink::numpify(const ICluster& cluster)
{
    // An empty cluster graph yields no arrays, hence no archive members,
    // and a downstream ClusterFileSource would then see EOS instead of an
    // empty cluster -- desynchronizing multi-stream consumers (e.g. the
    // live/dead pair feeding PointTreeBuilding).  The JSON serializer
    // always writes one (tiny) member per cluster and the source
    // dispatches per-member by extension, so fall back to it here.
    if (boost::num_vertices(cluster.graph()) == 0) {
        jsonify(cluster);
        return;
    }

    auto fn = [&](const std::string& name) {
        return this->fqn(cluster, name, "npy");
    };

    Aux::ClusterArrays::node_array_set_t nas;
    Aux::ClusterArrays::edge_array_set_t eas;
    Aux::ClusterArrays::to_arrays(cluster.graph(), nas, eas);

    // write nodes
    for (const auto& [nc,na] : nas) {
        std::string name = "_nodes";
        name[0] = nc;
        write_numpy(na, fn(name));
    }

    // write edges
    for (const auto& [ec, ea] : eas) {
        auto name = Aux::ClusterArrays::to_string(ec);
        name += "edges";
        write_numpy(ea, fn(name));
    }

}

void Sio::ClusterFileSink::configure(const WireCell::Configuration& cfg)
{
    m_outname = get(cfg, "outname", m_outname);
    m_out.clear();
    custard::output_filters(m_out, m_outname);
    if (m_out.empty()) {
        THROW(ValueError() << errmsg{"ClusterFileSink: unsupported outname: " + m_outname});
    }
    m_prefix = get<std::string>(cfg, "prefix", m_prefix);
    m_format = get<std::string>(cfg, "format", m_format);
    // force-override format if .npz
    if (String::endswith(m_outname, ".npz")) {
        m_format = "numpy";
    }
    if (m_format == "json") {
        m_serializer = [&](const ICluster& cluster){this->jsonify(cluster);};
    }
    else if (m_format == "dot") {
        m_serializer = [&](const ICluster& cluster){this->dotify(cluster);};
    }
    else if (m_format == "numpy") {
        m_serializer = [&](const ICluster& cluster){this->numpify(cluster);};
    }
    else if (m_format == "dummy") {
        m_serializer = [&](const ICluster& cluster){};
    }
    else {
        THROW(ValueError() << errmsg{"ClusterFileSink: unsupported format: " + m_format});
    }
}


bool Sio::ClusterFileSink::operator()(const ICluster::pointer& cluster)
{
    if (!cluster) {             // EOS
        log->debug("see EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    const auto& gr = cluster->graph();
    {
        // fixme: debugging.
        for (auto vtx : mir(boost::vertices(gr))) {
            const auto& node = gr[vtx];
            if (node.code() == 'b') {
                auto iblob = std::get<IBlob::pointer>(gr[vtx].ptr);
                Aux::BlobCategory bcat(iblob);
                if (bcat.ok()) continue;
                log->warn("malformed blob: {}", bcat.str());
            }
        }
    }
    log->debug("save cluster {} at call={}: {}", cluster->ident(), m_count, dumps(gr));

    m_serializer(*cluster);

    ++m_count;
    return true;
}
