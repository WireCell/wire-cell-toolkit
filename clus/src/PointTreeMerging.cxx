#include "WireCellClus/PointTreeMerging.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ExecMon.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"

WIRECELL_FACTORY(PointTreeMerging, WireCell::Clus::PointTreeMerging,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;


Clus::PointTreeMerging::PointTreeMerging()
    : Aux::Logger("PointTreeMerging", "clus")
{
}

std::vector<std::string> Clus::PointTreeMerging::input_types()
{
    log->debug("m_multiplicity {}", m_multiplicity);
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    log->debug("input_types: ret.size() {}", ret.size());
    return ret;
}

void Clus::PointTreeMerging::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);
    m_multiplicity = get<int>(cfg, "multiplicity", m_multiplicity);
    log->debug("{}", cfg);
    log->debug("m_multiplicity {}", m_multiplicity);
}

WireCell::Configuration Clus::PointTreeMerging::default_configuration() const
{
    Configuration cfg;
    return cfg;
}

void Clus::PointTreeMerging::finalize()
{
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
        log->debug("EOS at call {}", m_count++);
        return true;
    }
    if (neos) { raise<ValueError>("missing %d input tensors ", neos); }

    const int ident = invec[0]->ident();




    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) {
        inpath = String::format(inpath, ident);
    }
    auto root_live = as_pctree(*invec[0]->tensors(), inpath + "/live");
    if (!root_live) {
        log->error("Failed to get point cloud tree from \"{}\"", inpath);
        return false;
    }
    auto root_dead = as_pctree(*invec[0]->tensors(), inpath + "/dead");

    std::string outpath = m_outpath;
    if (outpath.find("%") != std::string::npos) {
        outpath = String::format(outpath, ident);
    }
    auto outtens = as_tensors(*root_live.get(), outpath + "/live");
    auto outtens_dead = as_tensors(*root_dead.get(), outpath + "/dead");
    for(const auto& ten : outtens) {
        log->debug("tensor {} {}", ten->metadata()["datapath"].asString(), ten->size());
        break;
    }
    outtens.insert(outtens.end(), outtens_dead.begin(), outtens_dead.end());
    outts = as_tensorset(outtens, ident);

    root_live = nullptr;
    root_dead = nullptr;

    m_count++;
    return true;
}