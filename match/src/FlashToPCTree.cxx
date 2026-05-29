#include "WireCellMatch/FlashToPCTree.h"

#include "WireCellSio/TensorFileSource.h"

#include "WireCellAux/TensorDMcommon.h"   // as_tensorset
#include "WireCellAux/TensorDMpointtree.h" // as_pctree, as_tensors

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointCloudArray.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/String.h"

WIRECELL_FACTORY(FlashToPCTree,
                 WireCell::Match::FlashToPCTree,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Match;

FlashToPCTree::FlashToPCTree() : Aux::Logger("FlashToPCTree", "match") {}
FlashToPCTree::~FlashToPCTree() = default;

void FlashToPCTree::configure(const WireCell::Configuration& cfg)
{
    m_input  = get(cfg, "input",  m_input);
    m_prefix = get(cfg, "prefix", m_prefix);
    m_inpath = get(cfg, "inpath", m_inpath);
    m_pcname = get(cfg, "pcname", m_pcname);

    if (m_input.empty()) {
        raise<ValueError>("FlashToPCTree: 'input' (flash tensor file) is required");
    }

    // Compose a TensorFileSource to read the raw flash tensor archive. Reusing
    // it keeps the parsed flash values byte-identical to the prior direct input.
    auto src = std::make_shared<Sio::TensorFileSource>();
    Configuration scfg = src->default_configuration();
    scfg["inname"] = m_input;
    scfg["prefix"] = m_prefix;
    src->configure(scfg);
    m_flashsrc = src;

    log->debug("FlashToPCTree: input={} prefix={} inpath={} pcname={}",
               m_input, m_prefix, m_inpath, m_pcname);
}

WireCell::Configuration FlashToPCTree::default_configuration() const
{
    Configuration cfg;
    cfg["input"]  = m_input;
    cfg["prefix"] = m_prefix;
    cfg["inpath"] = m_inpath;
    cfg["pcname"] = m_pcname;
    return cfg;
}

bool FlashToPCTree::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    const int ident = in->ident();
    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) inpath = String::format(inpath, ident);

    // Pull this event's flash tensor set from the composed source (one per
    // incoming pctree; the flash archive and cluster stream are event-aligned).
    ITensorSet::pointer flash_ts = nullptr;
    (*m_flashsrc)(flash_ts);

    // Deserialize the live tree and attach the flash measurement matrix as a
    // named point cloud on its root node.
    const auto& in_tens = *in->tensors();
    auto root_live = Aux::TensorDM::as_pctree(in_tens, inpath + "/live");
    if (!root_live) {
        raise<ValueError>("FlashToPCTree: no live pctree at '%s'", inpath);
    }

    PointCloud::Dataset ds;
    if (flash_ts && flash_ts->tensors() && !flash_ts->tensors()->empty()) {
        const auto& ten = flash_ts->tensors()->at(0);
        if (ten->shape().size() != 2) {
            raise<ValueError>("FlashToPCTree: flash tensor dim %d != 2", ten->shape().size());
        }
        const size_t nrow = ten->shape()[0];
        const size_t ncol = ten->shape()[1];
        // Store the [nflash, 1+nchan] matrix verbatim (col 0 = time,
        // cols 1.. = PE). The matcher rebuilds the same Opflash objects.
        PointCloud::Array arr((const double*) ten->data(),
                              PointCloud::Array::shape_t{nrow, ncol}, false);
        ds.add("value", std::move(arr));
        log->debug("FlashToPCTree: event {} attach flash [{},{}]", ident, nrow, ncol);
    }
    else {
        log->warn("FlashToPCTree: event {} has no flash tensor; attaching empty flash PC", ident);
    }
    root_live->value.local_pcs()[m_pcname] = std::move(ds);

    // Re-serialize live (now carrying the flash PC) and pass /dead through
    // unchanged, mirroring QLMatching's own output serialization.
    ITensor::vector outtens;
    auto tens_live = Aux::TensorDM::as_tensors(*root_live, inpath + "/live");
    outtens.insert(outtens.end(), tens_live.begin(), tens_live.end());

    auto root_dead = Aux::TensorDM::as_pctree(in_tens, inpath + "/dead");
    auto tens_dead = Aux::TensorDM::as_tensors(*root_dead, inpath + "/dead");
    outtens.insert(outtens.end(), tens_dead.begin(), tens_dead.end());

    out = Aux::TensorDM::as_tensorset(outtens, ident);
    ++m_count;
    return true;
}
