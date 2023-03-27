#include "WireCellImg/BlobDeclustering.h"

#include "WireCellAux/SimpleBlob.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/GraphTools.h"


WIRECELL_FACTORY(BlobDeclustering, WireCell::Img::BlobDeclustering,
                 WireCell::INamed,
                 WireCell::IBlobDeclustering)

using namespace WireCell;
using namespace WireCell::GraphTools;
using namespace WireCell::Img;
using namespace WireCell::Aux;

BlobDeclustering::BlobDeclustering()
    : Aux::Logger("BlobDeclustering", "img")
{
}

BlobDeclustering::~BlobDeclustering()
{
}

bool BlobDeclustering::operator()(const input_pointer& cluster, output_pointer& blobset)
{
    blobset = nullptr;
    if (!cluster) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    IBlob::vector blobs;
    const auto& gr = cluster->graph();
    for (const auto& vdesc : mir(boost::vertices(gr))) {
        const char code = gr[vdesc].code();
        if (code == 'b') {
            auto iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
            blobs.push_back(iblob);
        }
    }

    blobset = std::make_shared<SimpleBlobSet>(cluster->ident(), nullptr, blobs);

    log->debug("declustered {} blobs at call {}",
               blobs.size(), m_count++);
    return true;
}

        

