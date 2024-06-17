#include "WireCellClus/BlobSampling.h"
#include "WireCellUtil/PointTree.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"

WIRECELL_FACTORY(BlobSampling, WireCell::Clus::BlobSampling,
                 WireCell::INamed,
                 WireCell::IBlobSampling,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Tree;

BlobSampling::BlobSampling()
    : Aux::Logger("BlobSampling", "clus")
{
}


BlobSampling::~BlobSampling()
{
}


void BlobSampling::configure(const WireCell::Configuration& cfg)
{
    m_datapath = get(cfg, "datapath", m_datapath);
    auto samplers = cfg["samplers"];
    if (samplers.isNull()) {
        raise<ValueError>("add at least one entry to the \"samplers\" configuration parameter");
    }

    for (auto name : samplers.getMemberNames()) {
        auto tn = samplers[name].asString();
        if (tn.empty()) {
            raise<ValueError>("empty type/name for sampler \"%s\"", name);
        }
        log->debug("point cloud \"{}\" will be made by sampler \"{}\"",
                   name, tn);
        m_samplers[name] = Factory::find_tn<IBlobSampler>(tn); 
    }

}


WireCell::Configuration BlobSampling::default_configuration() const
{
    Configuration cfg;
    // eg:
    //    cfg["samplers"]["samples"] = "BlobSampler";
    cfg["datapath"] = m_datapath;
    return cfg;
}

bool BlobSampling::operator()(const input_pointer& blobset, output_pointer& tensorset)
{
    tensorset = nullptr;
    if (!blobset) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    auto iblobs = blobset->blobs();
    size_t nblobs = iblobs.size();

    Points::node_ptr root = std::make_unique<Points::node_t>();
    for (size_t bind=0; bind<nblobs; ++bind) {
        auto iblob = iblobs[bind];
        if (!iblob) {
            log->debug("skipping null blob {} of {}", bind, nblobs);
            continue;
        }
        named_pointclouds_t pcs;
        for (auto& [name, sampler] : m_samplers) {
            pcs.emplace(name, sampler->sample_blob(iblob, bind));
        }
        root->insert(Points(pcs));
    }
    // auto ds = m_sampler->sample_blobs(blobs);

    const int ident = blobset->ident();
    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }
    auto tens = as_tensors(*root.get(), datapath);
    tensorset = as_tensorset(tens, ident);

    log->debug("sampled {} blobs from set {} making {} tensors at call {}",
               nblobs, ident, tens.size(), m_count++);

    return true;
}

        

