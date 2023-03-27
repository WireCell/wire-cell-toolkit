#include "WireCellImg/BlobSetFanout.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include <iostream>

WIRECELL_FACTORY(BlobSetFanout, WireCell::Img::BlobSetFanout,
                 WireCell::INamed,
                 WireCell::IBlobSetFanout, WireCell::IConfigurable)

using namespace WireCell;
using namespace std;

Img::BlobSetFanout::BlobSetFanout(size_t multiplicity)
    : Aux::Logger("BlobSetFanout", "glue")
    , m_multiplicity(multiplicity)
{
}

Img::BlobSetFanout::~BlobSetFanout() {}

WireCell::Configuration Img::BlobSetFanout::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = (int) m_multiplicity;
    return cfg;
}

void Img::BlobSetFanout::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m <= 0) {
        log->critical("multiplicity must be positive, got {}", m);
        THROW(ValueError() << errmsg{"BlobSetFanout multiplicity must be positive"});
    }
    m_multiplicity = m;
}

std::vector<std::string> Img::BlobSetFanout::output_types()
{
    const std::string tname = std::string(typeid(output_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool Img::BlobSetFanout::operator()(const input_pointer& in, output_vector& outv)
{
    // Note: if "in" indicates EOS, just pass it on
    if (in) {
        const auto nblobs = in->blobs().size();
        // Too noisy for debug().
        log->debug("call={} fanout blob set {} with {}",
                   m_count, in->ident(), nblobs);
    }
    else {
        log->debug("EOS at call={}", m_count);
    }

    outv.resize(m_multiplicity, nullptr);
    for (size_t ind = 0; ind < m_multiplicity; ++ind) {
        outv[ind] = in;
    }

    ++m_count;
    return true;
}
