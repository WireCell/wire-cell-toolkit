#include "WireCellSio/ClusterFileSource.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/ClusterArrays.h" // to_cluster
// debugging
#include "WireCellUtil/GraphTools.h"
#include "WireCellAux/BlobTools.h"

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
using namespace WireCell::GraphTools;
using namespace WireCell::Sio;
using namespace WireCell::String;

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
    enum Type { bad, full, node, edge };
    enum Form { unknown, npy, json };

    Type type{bad};
    Form form{unknown};
    int ident{-1};
    std::string code{""};
};

// This throws ValueError if parse fails.  Catch in operator()
static
ParsedFilename parse_fname(WireCell::Log::logptr_t& log,
                           const std::string& fname,
                           const std::string& prefix)
{
    ParsedFilename ret;
    if (! startswith(fname, prefix)) {
        log->warn("file name {} does not start with prefix {}", fname, prefix);
        THROW(ValueError() << errmsg{"file name lacks prefix: " + fname});
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
    else {
        log->warn("file name {} has unknown extension", fname);
        THROW(ValueError() << errmsg{"unknown file name extension: " + fname});
        return ret;
    }

    auto parts = split(basename, "_");

    if (parts.size() < 2 or parts[0] != prefix) {
        log->warn("file name {} malformed", fname);
        THROW(ValueError() << errmsg{"malformed file name: " + fname});
        return ret;
    }

    const size_t nparts = parts.size();

    // cluster graph schema.
    // <prefix>_<clusterID>.json
    // or
    // <prefix>_<clusterID>_graph.json
    if (ret.form == ParsedFilename::json and (
            (nparts == 2 or (nparts == 3 and parts[2] == "graph")))) {
        ret.type = ParsedFilename::full;
        ret.ident = atoi(parts[1].c_str());
        // log->debug("file {} is cluster graph schema JSON", fname);
        return ret;
    }

    // cluster array schema.
    // <prefix>_<clusterID>_<nodecode>nodes.npy
    // <prefix>_<clusterID>_<edgecode>edges.npy
    if (nparts == 3 and ret.form == ParsedFilename::npy) {
        ret.ident = atoi(parts[1].c_str());
        const std::string aname = parts[2];
        if (aname.size() == 6 and endswith(aname, "nodes")) {
            ret.type = ParsedFilename::node;
            ret.code = aname.substr(0,1);
            // log->debug("file {} is cluster array schema Numpy of {} nodes", fname, ret.code);
            return ret;
        }
        if (aname.size() == 7 and endswith(aname, "edges")) {
            ret.type = ParsedFilename::edge;
            ret.code = aname.substr(0,2);
            // log->debug("file {} is cluster array schema Numpy of {} edges", fname, ret.code);
            return ret;
        }
    }
        
    log->warn("file {} is not supported", fname);
    THROW(ValueError() << errmsg{"file not supported: " + fname});
    // For cluster tensor schema: see TensorFileSource.
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
    cluster_graph_t cgraph = m_json_loader->load(jobj); // can throw
    return std::make_shared<SimpleCluster>(cgraph, ident);
}

ICluster::pointer ClusterFileSource::load_numpy(int ident)
{
    using namespace WireCell::Aux::ClusterArrays;
    node_array_set_t nas;
    edge_array_set_t eas;

    // Slurp many arrays.
    while (true) {
        if (! m_cur.fsize) {    // first pass filename is preloaded
            bool ok = load_filename();
            if (!ok) {
                break;
                // return nullptr;
            }                    
        }

        auto pf = parse_fname(log, m_cur.fname, m_prefix);
        if (ident != pf.ident) {
            break;              // next one is not ours
        }

        // load actual array
        pigenc::File pig;
        pig.read(m_in);
        auto shape = pig.header().shape();
        if (shape.size() != 2) {
            THROW(ValueError() << errmsg{"illegal shape"});
        }
        log->debug("file {} with type={} code={} ident={} shape=({},{})",
                   m_cur.fname, (int)pf.type, pf.code, ident, shape[0], shape[1]);
        
        if (pf.type == ParsedFilename::node) {
            const node_element_t* data = pig.as_type<node_element_t>();
            if (!data) {
                THROW(ValueError() << errmsg{"node array type mismatch"});
            }
            nas.emplace(pf.code[0], boost::const_multi_array_ref<node_element_t, 2>(data, shape));
        }
        else if (pf.type == ParsedFilename::edge) {
            const edge_element_t* data = pig.as_type<edge_element_t>();
            if (!data) {
                THROW(ValueError() << errmsg{"edge array type mismatch"});
            }
            auto ecode = edge_code(pf.code[0], pf.code[1]);
            eas.emplace(ecode, boost::const_multi_array_ref<edge_element_t, 2>(data, shape));
        }
        else {
            this->clear_load();
        }
        clear();
            
        /// originally, we had the number of arrays hard wired to be a full
        /// ICluster.  When running, eg, BlobClustering but BlobGrouping there
        /// may not be any "m" nodes.  But, in general, make judgment, best to
        /// delay judgement.
        // if (nas.size() == 5 and eas.size() == 7) {
        //     // log->debug("completed numpy load for ident={}", ident);
        //     break;
        // }
    };

    //// See above comments
    // if (nas.size() != 5) {
    //     log->error("ident={} failed to load all node arrays, loaded {}:", ident, nas.size());
    //     for (const auto& [code,arr] : nas) {
    //         log->error("\t{}: ({},{})", code, arr.shape()[0], arr.shape()[1]);
    //     }
    //     return nullptr;
    // }
    // if (eas.size() != 7) {
    //     log->error("ident={} failed to load all edge arrays, loaded {}:", ident, eas.size());
    //     for (const auto& [code,arr] : eas) {
    //         log->error("\t{}: ({},{})", code, arr.shape()[0], arr.shape()[1]);
    //     }
    //     return nullptr;
    // }

    auto graph = to_cluster(nas, eas, m_anodes);
    return std::make_shared<SimpleCluster>(graph, ident);
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

    for (const auto& tn : cfg["anodes"]) {
        auto ianode = Factory::find_tn<IAnodePlane>(tn.asString());
        m_anodes.push_back(ianode);
    }
    if (m_anodes.empty()) {
        THROW(ValueError() << errmsg{"ClusterFileSource: no anodes given"});
    }
    m_json_loader = std::make_unique<ClusterLoader>(m_anodes);
}


bool ClusterFileSource::load_filename()
{
    custard::read(m_in, m_cur.fname, m_cur.fsize);
    if (m_in.eof()) {
        // log->debug("eof in file: {} size {}", m_cur.fname, m_cur.fsize);
        return false;
    }
    if (!m_in) {
        log->critical("call={}, read stream error with file={}",
                      m_count, m_inname);
        return false;
    }
    if (!m_cur.fsize) {
        log->critical("call={}, short read from file={}",
                      m_count, m_inname);
        return false;
    }

    // log->debug("loading {} of size {}", m_cur.fname, m_cur.fsize);

    return true;
}

ICluster::pointer ClusterFileSource::dispatch()
{
    auto pf = parse_fname(log, m_cur.fname, m_prefix);
    // log->debug("dispatch file {} with type={} code={} ident={}",
    //            m_cur.fname, pf.type, pf.code, pf.ident);
    if (pf.form == ParsedFilename::json) {
        return load_json(pf.ident);
    }
    if (endswith(m_cur.fname, ".npy")) {
        return load_numpy(pf.ident);
    }
    log->warn("do not know how to dispatch file {} with type={} code={} ident={}",
              m_cur.fname, (int)pf.type, pf.code, pf.ident);
    return nullptr;
}

void ClusterFileSource::clear_load()
{
    log->warn("call={} skipping unsupported file \"{}\" of size {} in stream: {}",
              m_count, m_cur.fname, m_cur.fsize, m_inname);
    m_in.seekg(m_cur.fsize, m_in.cur);
    clear();
}

ICluster::pointer ClusterFileSource::load()
{
    while (true) {
        if (! load_filename()) {
            return nullptr;
        }
        auto ret = dispatch();
        if (ret) return ret;

        this->clear_load();
    }
}

bool ClusterFileSource::operator()(ICluster::pointer& cluster)
{
    cluster = nullptr;
    if (m_eos_sent) {
        return false;
    }
    try {
        cluster = load();
    }
    catch (ValueError& err) {
        log->error(err.what());
        log->error("failed to load file, ending stream at call={}", m_count);
        cluster = nullptr;
        m_eos_sent = true;
        return false;
    }
        
    if (! cluster) {
        log->debug("see EOS at call={}", m_count);
        m_eos_sent = true;
    }
    else {
        // fixme: debugging. 
        const auto& cgraph = cluster->graph();
        for (auto vtx : mir(boost::vertices(cgraph))) {
            const auto& node = cgraph[vtx];
            if (node.code() == 'b') {
                auto iblob = std::get<IBlob::pointer>(cgraph[vtx].ptr);
                Aux::BlobCategory bcat(iblob);
                if (bcat.ok()) continue;
                log->warn("malformed blob: {}", bcat.str());
            }
        }

        log->debug("load cluster {} at call={}: {}", cluster->ident(), m_count, dumps(cluster->graph()));
    }
    ++m_count;
    return true;
}

void ClusterFileSource::clear()
{
    m_cur = header_t();
}
