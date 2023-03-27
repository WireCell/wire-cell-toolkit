#include "WireCellSio/ClusterFileSink.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/ClusterArrays.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/custard/custard_boost.hpp"

WIRECELL_FACTORY(ClusterFileSink, WireCell::Sio::ClusterFileSink,
                 WireCell::INamed,
                 WireCell::IClusterSink, WireCell::ITerminal,
                 WireCell::IConfigurable)

static void dump_cluster(WireCell::Log::logptr_t& log,
                         const WireCell::ICluster::pointer& cluster)
{
    const auto& gr = cluster->graph();
    const auto known = WireCell::cluster_node_t::known_codes;
    const size_t ncodes = known.size();
    std::vector<size_t> counts(ncodes, 0);

    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& vobj = gr[vtx];
        ++counts[vobj.ptr.index() - 1]; // node type index=0 is the bogus node type
    }

    std::stringstream ss;
    for (size_t ind=0; ind < ncodes; ++ind) {
        ss << known[ind] << ":" << counts[ind] << " ";
    }

    log->debug("cluster={} nvertices={} nedges={} types: {}",
               cluster->ident(),
               boost::num_vertices(gr),
               boost::num_edges(gr),
               ss.str());
}

using namespace WireCell;

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
    auto fn = [&](const std::string& name) {
        return this->fqn(cluster, name, "npy");
    };

    Aux::ClusterArrays ca(cluster.graph());

    // write nodes
    for (auto nc : ca.node_codes()) {
        const auto& na = ca.node_array(nc);
        std::string name = "_nodes";
        name[0] = nc;
        write_numpy(na, fn(name));
    }

    // write edges
    for (auto ec : ca.edge_codes()) {
        const auto& ea = ca.edge_array(ec);
        auto name = ca.edge_code_str(ec);
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

    dump_cluster(log, cluster);
    m_serializer(*cluster);

    ++m_count;
    return true;
}
