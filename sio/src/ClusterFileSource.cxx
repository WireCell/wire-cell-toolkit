#include "WireCellSio/ClusterFileSource.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/ClusterArrays.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleCluster.h"

#include "WireCellUtil/custard/custard_boost.hpp"

WIRECELL_FACTORY(ClusterFileSource, WireCell::Sio::ClusterFileSource,
                 WireCell::INamed,
                 WireCell::IClusterSource, WireCell::ITerminal,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Sio;
using namespace WireCell::String;

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

ClusterFileSource::ClusterFileSource()
    : Logger("ClusterFileSource", "io")
{
}

ClusterFileSource::~ClusterFileSource()
{
}

void ClusterFileSource::finalize()
{
}

WireCell::Configuration ClusterFileSource::default_configuration() const
{
    Configuration cfg;

    cfg["inname"] = m_inname;
    cfg["prefix"] = m_prefix;
    cfg["anodes"] = Json::arrayValue;

    return cfg;
}

struct ParsedFilename {
    enum Type { bad, full, node, edge, attr };
    enum Form { unknown, npy, json };
    Type type{bad};
    Form form{unknown};
    int ident{-1};
    std::string atype{""}, aname{""};
};
static
ParsedFilename parse_fname(const std::string& fname,
                           const std::string& prefix)
{
    ParsedFilename ret;
    if (! startswith(fname, prefix)) {
        return ret;
    }

    std::string basename;
    if (endswith(fname, ".json")) {
        ret.form = ParsedFilename::json;
        basename = fname.substr(0, fname.size() - 5);
    }
    else if (endswith(fname, ".npy")) {
        ret.form = ParsedFilename::npy;
        basename = fname.substr(0, fname.size() - 4);
    }
    else return ret;

    auto parts = split(basename, "_");

    if (parts.size() < 2 or parts[0] != prefix) {
        return ret;
    }

    //<prefix>_<clusterID>.json
    if (parts.size() == 2 and ret.form == ParsedFilename::json) {
        ret.type = ParsedFilename::full;
        ret.ident = atoi(parts[1].c_str());
        return ret;
    }

    // <prefix>_<clusterID>_node_<type>.npy
    // <prefix>_<clusterID>_edge_<type>.npy
    // <prefix>_<clusterID>_node_<type>_<attrname>.npy
    // <prefix>_<clusterID>_edge_<type>_<attrname>.npy
    if (parts.size() > 3 and ret.form == ParsedFilename::npy) {
        if (parts[2] == "node") {
            ret.type = ParsedFilename::node;
        }
        else if (parts[2] == "edge") {
            ret.type = ParsedFilename::edge;
        }
        else {
        }
        if (parts.size() == 5) {
            ret.type = ParsedFilename::attr;
            ret.atype = parts[3];
            ret.aname = parts[4];
        }
    }
    return ret;
}

ICluster::pointer ClusterFileSource::load_json(int ident)
{
    Configuration jobj;

    // Read into a buffer string to assure we consume fsize bytes and
    // protect against malformed JSON or leading/trailing spaces.
    std::string buffer;
    buffer.resize(m_cur.fsize);
    m_in.read(buffer.data(), buffer.size());
    clear();
    if (!m_in) {
        return nullptr;
    }
    std::istringstream ss(buffer);
    ss >> jobj;
    auto cgraph = m_json_loader->load(jobj);
    return std::make_shared<SimpleCluster>(cgraph, ident);
}

ICluster::pointer ClusterFileSource::load_numpy(int ident)
{
    log->warn("no numpy support yet, skipping {} at call={}",
              m_cur.fname, m_count);
    return nullptr;
}

void ClusterFileSource::configure(const WireCell::Configuration& cfg)
{
    m_inname = get(cfg, "inname", m_inname);
    m_in.clear();
    custard::input_filters(m_in, m_inname);
    if (m_in.empty()) {
        THROW(ValueError() << errmsg{"ClusterFileSource: unsupported inname: " + m_inname});
    }
    m_prefix = get<std::string>(cfg, "prefix", m_prefix);
    std::vector<IAnodePlane::pointer> anodes;
    for (const auto& tn : cfg["anodes"]) {
        auto ianode = Factory::find_tn<IAnodePlane>(tn.asString());
        anodes.push_back(ianode);
    }
    if (anodes.empty()) {
        THROW(ValueError() << errmsg{"ClusterFileSource: no anodes given"});
    }
    m_json_loader = std::make_unique<ClusterLoader>(anodes);
}


ICluster::pointer ClusterFileSource::load()
{
    while (true) {
        custard::read(m_in, m_cur.fname, m_cur.fsize);
        if (m_in.eof()) {
            return nullptr;
        }
        if (!m_in) {
            log->critical("call={}, read stream error with file={}",
                          m_count, m_inname);
            return nullptr;
        }
        if (!m_cur.fsize) {
            log->critical("call={}, short read from file={}",
                          m_count, m_inname);
            return nullptr;
        }

        auto pf = parse_fname(m_cur.fname, m_prefix);
        if (pf.form == ParsedFilename::json) {
            return load_json(pf.ident);
        }
        if (endswith(m_cur.fname, ".npy")) {
            return load_numpy(pf.ident);
        }
        log->warn("call={} skipping unsupported file {} in stream {}",
                  m_count, m_cur.fname, m_inname);
        m_in.seekg(m_cur.fsize, m_in.cur);
        clear();
    }
}

bool ClusterFileSource::operator()(ICluster::pointer& cluster)
{
    cluster = nullptr;
    if (m_eos_sent) {
        return false;
    }
    cluster = load();
    if (! cluster) {
        log->debug("see EOS at call={}", m_count);
        m_eos_sent = true;
    }
    else {
        log->debug("load cluster at call={}", m_count);
    }
    ++m_count;
    return true;
}

void ClusterFileSource::clear()
{
    m_cur = header_t();
}
