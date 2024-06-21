#include "WireCellImg/BlobSetMerge.h"
#include "WireCellImg/ImgData.h"
#include "WireCellAux/SimpleBlob.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(BlobSetMerge, WireCell::Img::BlobSetMerge,
                 WireCell::INamed,
                 WireCell::IBlobSetFanin, WireCell::IConfigurable)
using namespace WireCell;

Img::BlobSetMerge::BlobSetMerge()
    : Aux::Logger("BlobSetMerge", "glue")
    , m_multiplicity(0)
{
}

Img::BlobSetMerge::~BlobSetMerge() {}

WireCell::Configuration Img::BlobSetMerge::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = (int) m_multiplicity;
    return cfg;
}

void Img::BlobSetMerge::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m <= 0) {
        THROW(ValueError() << errmsg{"BlobSetMerge multiplicity must be positive"});
    }
    m_multiplicity = m;
}

std::vector<std::string> Img::BlobSetMerge::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool Img::BlobSetMerge::operator()(const input_vector& invec, output_pointer& out)
{
    Img::Data::Slice *oslice = nullptr;
    auto sbs = std::make_shared<Aux::SimpleBlobSet>(0, nullptr);

    // TODO assign new blob ident
    // reset the counter with EOS

    int neos = 0;
    std::stringstream perset;
    for (const auto& ibs : invec) {
        if (!ibs) {
            ++neos;
            break;
        }
        SPDLOG_LOGGER_TRACE(log,"in: slice ident: {} time: {} activites: {} blobs: {}",
                            ibs->slice()->ident(), ibs->slice()->start(),
                            ibs->slice()->activity().size(), ibs->blobs().size());
        ISlice::pointer newslice = ibs->slice();
        if (!sbs->slice()) {
            oslice = new Img::Data::Slice(newslice->frame(), newslice->ident(),
                                          newslice->start(), newslice->span());
            sbs->m_slice = ISlice::pointer(oslice);
            sbs->m_ident = newslice->ident();
        }
        // make sure times are sync'ed
        if(sbs->slice()->start()!=newslice->start()) {
            THROW(RuntimeError() << errmsg{"blobs not in sync'ed slices"});
        }
        oslice->merge(newslice);
        for (const auto& iblob : ibs->blobs()) {
            // FIXME: ident can be the same?
            Aux::SimpleBlob* sb = new Aux::SimpleBlob(iblob->ident(), iblob->value(),
                                                      iblob->uncertainty(),
            iblob->shape(), sbs->m_slice, iblob->face());
            sbs->m_blobs.push_back(IBlob::pointer(sb));
        }
        perset << " " << ibs->blobs().size();

        SPDLOG_LOGGER_TRACE(log, "out: slice ident: {} time: {} activites: {} blobs: {}",
                            sbs->slice()->ident(), sbs->slice()->start(),
                            sbs->slice()->activity().size(), sbs->blobs().size());
    }
    perset << " ";
    if (neos) {
        out = nullptr;
        log->debug("EOS");
        return true;
    }

    out = sbs;

    log->trace("sync'ed {} blobs: {}", sbs->m_blobs.size(), perset.str());
    return true;
}
