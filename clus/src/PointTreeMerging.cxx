#include "WireCellClus/PointTreeMerging.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"

#include <map>

WIRECELL_FACTORY(PointTreeMerging, WireCell::Clus::PointTreeMerging,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;


Clus::PointTreeMerging::PointTreeMerging()
    : Aux::Logger("PointTreeMerging", "clus")
{
}

std::vector<std::string> Clus::PointTreeMerging::input_types()
{
    SPDLOG_LOGGER_TRACE(log, "m_multiplicity {}", m_multiplicity);
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    SPDLOG_LOGGER_TRACE(log, "input_types: ret.size() {}", ret.size());
    return ret;
}

void Clus::PointTreeMerging::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);
    m_multiplicity = get<int>(cfg, "multiplicity", m_multiplicity);
    m_tolerate_missing = get<bool>(cfg, "tolerate_missing", m_tolerate_missing);
    SPDLOG_LOGGER_TRACE(log, "{}", cfg);
    SPDLOG_LOGGER_TRACE(log, "m_multiplicity {}", m_multiplicity);
}

WireCell::Configuration Clus::PointTreeMerging::default_configuration() const
{
    Configuration cfg;
    return cfg;
}

void Clus::PointTreeMerging::finalize()
{
}

static size_t normalize_pctree_local_pcs(Points::node_t* root)
{
    if (!root) {
        return 0;
    }

    std::vector<Points::node_t*> nodes;
    for (auto& noderef : root->depth()) {
        nodes.push_back(&noderef);
    }

    std::map<std::string, std::map<std::string, WireCell::PointCloud::Array>> templates;
    for (const auto* node : nodes) {
        for (const auto& [pcname, ds] : node->value.local_pcs()) {
            auto& pc_templates = templates[pcname];
            for (const auto& key : ds.keys()) {
                if (pc_templates.count(key)) {
                    continue;
                }
                pc_templates.emplace(key, *ds.get(key));
            }
        }
    }

    size_t nadded = 0;
    for (auto* node : nodes) {
        for (auto& [pcname, ds] : node->value.local_pcs()) {
            const auto it = templates.find(pcname);
            if (it == templates.end()) {
                continue;
            }
            for (const auto& [key, arr_template] : it->second) {
                if (ds.has(key)) {
                    continue;
                }
                ds.add(key, arr_template.zeros_like(ds.size_major()));
                ++nadded;
            }
        }
    }
    return nadded;
}

static void merge_pct(Points::node_t* tgt, Points::node_t* src)
{
    if (!src) {
        return;
    }

    // merge local pcs for root node
    auto tgt_pc = tgt->value.local_pcs();
    for (const auto& src_pc : src->value.local_pcs()) {
        auto name = src_pc.first;
        if (tgt_pc.find(name) == tgt_pc.end()) {
            tgt_pc.emplace(name, src_pc.second);
        } else {
            auto& tgt_pcds = tgt_pc[name];
            tgt_pcds.append(src_pc.second);
        }
    }

    // merge children
    bool notify_value = true;
    tgt->take_children(*src, notify_value);
}


bool Clus::PointTreeMerging::operator()(const input_vector& invec, output_pointer& outts)
{
    outts = nullptr;
    if (invec.empty()) {
        raise<ValueError>("no input tensors");
        return true;
    }
    // check input size
    if (invec.size() != m_multiplicity) {
        raise<ValueError>("unexpected multiplicity got %d want %d", invec.size(), m_multiplicity);
        return true;
    }
    // boilerplate for EOS handling
    size_t neos = 0;
    for (const auto& in : invec) {
        if (!in) { ++neos; }
    }
    if (neos == invec.size()) {
        // all inputs are EOS, good.
        SPDLOG_LOGGER_TRACE(log, "EOS at call {}", m_count++);
        return true;
    }
    if (neos) { raise<ValueError>("missing %d input tensors ", neos); }

    const int ident = invec[0]->ident();

    // Format m_inpath with a given ident (if m_inpath contains "%").
    auto format_inpath = [&](int id) {
        std::string p = m_inpath;
        if (p.find("%") != std::string::npos) {
            p = String::format(p, id);
        }
        return p;
    };

    // Per-input path: when m_tolerate_missing, use each input's own ident so
    // tensor sets emitted with different idents (e.g. empty placeholder sets
    // carry ident=0 while real sets carry the art event number) are still
    // found at their actual datapath. Otherwise keep legacy behavior of using
    // invec[0]'s ident for all inputs.
    auto inpath_for = [&](const ITensorSet::pointer& ts) {
        return format_inpath(m_tolerate_missing ? ts->ident() : ident);
    };

    // Wrap as_pctree: when m_tolerate_missing is set, a missing datapath
    // (raised as KeyError by TensorIndex::at) is treated as an empty tree.
    auto as_pctree_tol = [&](const ITensor::vector& tens, const std::string& path)
        -> std::unique_ptr<WireCell::PointCloud::Tree::Points::node_t> {
        try {
            return as_pctree(tens, path);
        }
        catch (const KeyError& e) {
            if (!m_tolerate_missing) throw;
            SPDLOG_LOGGER_DEBUG(log, "tolerating missing datapath '{}' (empty tree)", path);
            return std::make_unique<WireCell::PointCloud::Tree::Points::node_t>();
        }
    };

    auto root_live = as_pctree_tol(*invec[0]->tensors(), inpath_for(invec[0]) + "/live");
    if (!root_live) {
        SPDLOG_LOGGER_ERROR(log, "Failed to get point cloud tree from \"{}\"", inpath_for(invec[0]) + "/live");
        return false;
    }
    auto root_dead = as_pctree_tol(*invec[0]->tensors(), inpath_for(invec[0]) + "/dead");
    if (!root_dead) {
        SPDLOG_LOGGER_ERROR(log, "Failed to get point cloud tree from \"{}\"", inpath_for(invec[0]) + "/dead");
        return false;
    }

    SPDLOG_LOGGER_DEBUG(log, "input[0] ident={} path='{}' live children={}",
                        invec[0]->ident(), inpath_for(invec[0]), root_live->nchildren());
    // merge
    for (size_t i = 1; i < invec.size(); ++i) {
        if (!invec[i]) {
            raise<ValueError>("missing input tensor %d", i);
        }
        auto src_live = as_pctree_tol(*invec[i]->tensors(), inpath_for(invec[i]) + "/live");
        auto src_dead = as_pctree_tol(*invec[i]->tensors(), inpath_for(invec[i]) + "/dead");
        SPDLOG_LOGGER_DEBUG(log, "input[{}] ident={} path='{}' live children={}",
                            i, invec[i]->ident(), inpath_for(invec[i]),
                            src_live ? src_live->nchildren() : 0);
        merge_pct(root_live.get(), src_live.get());
        merge_pct(root_dead.get(), src_dead.get());
    }


    SPDLOG_LOGGER_DEBUG(log, "merged live PC tree with {} children", root_live->nchildren());
    SPDLOG_LOGGER_DEBUG(log, "merged dead PC tree with {} children", root_dead->nchildren());
    const auto nlive_added = normalize_pctree_local_pcs(root_live.get());
    const auto ndead_added = normalize_pctree_local_pcs(root_dead.get());
    SPDLOG_LOGGER_DEBUG(log, "normalized merged PC trees with {} live and {} dead arrays added",
                        nlive_added, ndead_added);

    // output
    std::string outpath = m_outpath;
    if (outpath.find("%") != std::string::npos) {
        outpath = String::format(outpath, ident);
    }
    auto outtens = as_tensors(*root_live.get(), outpath + "/live");
    auto outtens_dead = as_tensors(*root_dead.get(), outpath + "/dead");
    outtens.insert(outtens.end(), outtens_dead.begin(), outtens_dead.end());
    for(const auto& ten : outtens) {
        SPDLOG_LOGGER_TRACE(log, "outtens {} {}", ten->metadata()["datapath"].asString(), ten->size());
        break;
    }
    outts = as_tensorset(outtens, ident);

    root_live = nullptr;
    root_dead = nullptr;

    m_count++;
    return true;
}
