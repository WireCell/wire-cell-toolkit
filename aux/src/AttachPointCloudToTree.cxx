#include "WireCellAux/AttachPointCloudToTree.h"

#include "WireCellAux/TensorDMcommon.h"    // as_tensorset
#include "WireCellAux/TensorDMpointtree.h" // as_pctree, as_tensors

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointCloudArray.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/String.h"

WIRECELL_FACTORY(AttachPointCloudToTree,
                 WireCell::Aux::AttachPointCloudToTree,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;

Aux::AttachPointCloudToTree::AttachPointCloudToTree()
    : Aux::Logger("AttachPointCloudToTree", "aux")
{
}

Aux::AttachPointCloudToTree::~AttachPointCloudToTree() = default;

void Aux::AttachPointCloudToTree::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_pcname = get(cfg, "pcname", m_pcname);

    if (m_pcname.empty()) {
        raise<ValueError>("AttachPointCloudToTree: 'pcname' (root-node PC name) is required");
    }

    log->debug("AttachPointCloudToTree: inpath={} pcname={}", m_inpath, m_pcname);
}

WireCell::Configuration Aux::AttachPointCloudToTree::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"] = m_inpath;
    cfg["pcname"] = m_pcname;
    return cfg;
}

std::vector<std::string> Aux::AttachPointCloudToTree::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

bool Aux::AttachPointCloudToTree::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;

    // Inputs are framework-synced: any EOS is EOS.
    for (const auto& in : invec) {
        if (!in) {
            log->debug("EOS at call {}", m_count++);
            return true;
        }
    }

    const auto& tree_ts = invec[0]; // pctree (live + dead)
    const auto& data_ts = invec[1]; // tensor set to attach

    const int ident = tree_ts->ident();
    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) inpath = String::format(inpath, ident);

    // Deserialize the live tree and attach the data tensor as a named PC on the
    // root node.
    const auto& tree_tens = *tree_ts->tensors();
    auto root_live = Aux::TensorDM::as_pctree(tree_tens, inpath + "/live");
    if (!root_live) {
        raise<ValueError>("AttachPointCloudToTree: no live pctree at '%s'", inpath);
    }

    PointCloud::Dataset ds;
    const auto& data_tens = data_ts->tensors();
    if (data_tens && !data_tens->empty()) {
        const auto& ten = data_tens->at(0);
        if (ten->shape().size() != 2) {
            raise<ValueError>("AttachPointCloudToTree: tensor dim %d != 2", ten->shape().size());
        }
        const size_t nrow = ten->shape()[0];
        const size_t ncol = ten->shape()[1];
        // Store the [nrow, ncol] matrix verbatim; the downstream consumer
        // interprets it.
        PointCloud::Array arr((const double*) ten->data(),
                              PointCloud::Array::shape_t{nrow, ncol}, false);
        ds.add("value", std::move(arr));
        log->debug("AttachPointCloudToTree: ident {} attach pc \"{}\" [{},{}]",
                   ident, m_pcname, nrow, ncol);
    }
    else {
        log->warn("AttachPointCloudToTree: ident {} no tensor on port 1; attaching empty pc \"{}\"",
                  ident, m_pcname);
    }
    root_live->value.local_pcs()[m_pcname] = std::move(ds);

    // Re-serialize live (now carrying the new PC) and pass /dead through.
    ITensor::vector outtens;
    auto tens_live = Aux::TensorDM::as_tensors(*root_live, inpath + "/live");
    outtens.insert(outtens.end(), tens_live.begin(), tens_live.end());

    auto root_dead = Aux::TensorDM::as_pctree(tree_tens, inpath + "/dead");
    auto tens_dead = Aux::TensorDM::as_tensors(*root_dead, inpath + "/dead");
    outtens.insert(outtens.end(), tens_dead.begin(), tens_dead.end());

    out = Aux::TensorDM::as_tensorset(outtens, ident);
    ++m_count;
    return true;
}
