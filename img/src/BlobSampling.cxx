#include "WireCellImg/BlobSampling.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/IndexedGraph.h"
#include "WireCellAux/TensorDM.h"

WIRECELL_FACTORY(BlobSampling, WireCell::Img::BlobSampling,
                 WireCell::INamed,
                 WireCell::IBlobSampling,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;

BlobSampling::BlobSampling()
    : Aux::Logger("BlobSampling", "img")
{
}


BlobSampling::~BlobSampling()
{
}


void BlobSampling::configure(const WireCell::Configuration& cfg)
{
    std::string tn = get<std::string>(cfg, "sampler", "BlobSampler");
    m_sampler = Factory::find_tn<IBlobSampler>(tn);
    m_datapath = get(cfg, "datapath", m_datapath);
}


WireCell::Configuration BlobSampling::default_configuration() const
{
    Configuration cfg;
    cfg["sampler"] = "BlobSampler";
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

    const int ident = blobset->ident();
    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }

    auto blobs = blobset->blobs();
    auto ds = m_sampler->sample_blobs(blobs);
    auto tens = as_tensors(ds, datapath);
    tensorset = as_tensorset(tens, ident);

    log->debug("sampled {} blobs from set {} making {} points at call {}",
               blobs.size(), ident, ds.size_major(), m_count++);
    return true;
}

        

