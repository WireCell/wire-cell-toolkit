#include "WireCellFlash/OpHitMerge.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include <algorithm>

WIRECELL_FACTORY(OpHitMerge, WireCell::Flash::OpHitMerge,
                 WireCell::INamed,
                 WireCell::IConfigurable,
                 WireCell::ITensorSetFanin)

using namespace WireCell;

Flash::OpHitMerge::OpHitMerge()
  : Aux::Logger("OpHitMerge", "flash")
{
}

Flash::OpHitMerge::~OpHitMerge() {}

WireCell::Configuration Flash::OpHitMerge::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = m_multiplicity;
    cfg["meta_port"] = m_meta_port;
    return cfg;
}

void Flash::OpHitMerge::configure(const WireCell::Configuration& cfg)
{
    m_multiplicity = get(cfg, "multiplicity", m_multiplicity);
    m_meta_port = get(cfg, "meta_port", m_meta_port);
    if (m_multiplicity < 1) {
        raise<ValueError>("OpHitMerge: illegal multiplicity %d", m_multiplicity);
    }
    if (m_meta_port < 0 || m_meta_port >= m_multiplicity) {
        raise<ValueError>("OpHitMerge: illegal meta_port %d", m_meta_port);
    }
}

std::vector<std::string> Flash::OpHitMerge::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

// Locate the "ophits" tensor in a set (by metadata name, else tensor 0) --
// mirrors OpFlashFinder::operator().
static ITensor::pointer find_ophits(const ITensorSet::pointer& ts)
{
    for (const auto& t : *ts->tensors()) {
        if (t->metadata()["name"].asString() == "ophits") {
            return t;
        }
    }
    return ts->tensors()->at(0);
}

bool Flash::OpHitMerge::operator()(const input_vector& in, output_pointer& out)
{
    out = nullptr;
    if (in.empty()) {
        return true;  // eos
    }
    for (const auto& ts : in) {
        if (!ts) {        // inputs are synced, so any EOS is EOS
            return true;
        }
    }

    // Gather the ophits tensors and total their rows.
    std::vector<ITensor::pointer> parts;
    size_t ncol = 0, nrow_total = 0;
    for (const auto& ts : in) {
        ITensor::pointer t = find_ophits(ts);
        const auto shape = t->shape();
        if (ncol == 0) {
            ncol = shape[1];
        }
        else if (shape[1] != ncol) {
            raise<ValueError>("OpHitMerge: ophits column mismatch %d vs %d",
                              (int) shape[1], (int) ncol);
        }
        nrow_total += shape[0];
        parts.push_back(t);
    }

    // Row-concatenate in port order.
    std::vector<double> merged(nrow_total * ncol);
    size_t off = 0;
    for (const auto& t : parts) {
        const size_t n = t->shape()[0] * ncol;
        const double* d = (const double*) t->data();
        std::copy(d, d + n, merged.data() + off);
        off += n;
    }

    ITensor::vector* tensors = new ITensor::vector;
    {
        Configuration md;
        md["name"] = "ophits";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nrow_total, ncol}, merged.data(), md));
    }

    const auto& meta_ts = in[m_meta_port];
    out = std::make_shared<Aux::SimpleTensorSet>(meta_ts->ident(), meta_ts->metadata(),
                                                 ITensor::shared_vector(tensors));
    log->debug("merged {} ophits tensors -> {} hits x {} cols",
               in.size(), nrow_total, ncol);
    return true;
}
