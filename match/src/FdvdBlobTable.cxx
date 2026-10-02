#include "WireCellMatch/FdvdBlobTable.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/NamedFactory.h"

#include <algorithm>
#include <numeric>

WIRECELL_FACTORY(FdvdBlobTable, WireCell::Match::FdvdBlobTable,
                 WireCell::INamed,
                 WireCell::IClusterFaninTensorSet, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::GraphTools::mir;

Match::FdvdBlobTable::FdvdBlobTable()
  : Aux::Logger("FdvdBlobTable", "match")
{
}

Match::FdvdBlobTable::~FdvdBlobTable() {}

std::vector<std::string> Match::FdvdBlobTable::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

WireCell::Configuration Match::FdvdBlobTable::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = (int) m_multiplicity;
    cfg["tick"] = m_tick;
    cfg["xw_mm"] = m_xw_mm;
    cfg["v_mm_per_us"] = m_v_mm_per_us;
    cfg["toff_us"] = m_toff_us;
    cfg["stored_corners"] = m_stored_corners;
    return cfg;
}

void Match::FdvdBlobTable::configure(const WireCell::Configuration& cfg)
{
    const int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m <= 0) raise<ValueError>("FdvdBlobTable: multiplicity must be positive");
    m_multiplicity = m;
    m_labels.clear();
    for (const auto& l : cfg["labels"]) m_labels.push_back(l.asString());
    if (!m_labels.empty() && m_labels.size() != m_multiplicity) {
        raise<ValueError>("FdvdBlobTable: %d labels for multiplicity %d", (int) m_labels.size(), (int) m_multiplicity);
    }
    m_tick = get(cfg, "tick", m_tick);
    m_xw_mm = get(cfg, "xw_mm", m_xw_mm);
    m_v_mm_per_us = get(cfg, "v_mm_per_us", m_v_mm_per_us);
    m_toff_us = get(cfg, "toff_us", m_toff_us);
    m_stored_corners = get(cfg, "stored_corners", m_stored_corners);
}

bool Match::FdvdBlobTable::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    size_t neos = 0;
    for (const auto& in : invec) neos += !in;
    if (neos == invec.size()) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    if (neos) raise<ValueError>("FdvdBlobTable: %d of %d inputs at EOS", (int) neos, (int) invec.size());
    if (invec.size() != m_multiplicity) {
        raise<ValueError>("FdvdBlobTable: %d inputs, multiplicity %d", (int) invec.size(), (int) m_multiplicity);
    }
    std::vector<size_t> ports(invec.size());
    std::iota(ports.begin(), ports.end(), 0);
    if (!m_labels.empty()) {
        std::stable_sort(ports.begin(), ports.end(), [&](size_t a, size_t b) { return m_labels[a] < m_labels[b]; });
    }
    std::vector<double> rows;
    size_t order = 0;
    m_nstored = 0;
    for (size_t port : ports) {
        const auto& gr = invec[port]->graph();
        for (const auto& vdesc : mir(boost::vertices(gr))) {
            if (gr[vdesc].code() != 'b') continue;
            const auto iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
            const auto face = iblob->face();
            const auto& coords = face->raygrid();
            const auto& shape = iblob->shape();
            const auto islice = iblob->slice();
            const double start = islice->start(), span = islice->span();
            double row[NCOL] = {0};
            row[0] = (double) order++;
            row[1] = face->anode();
            row[2] = face->which();
            row[3] = (int) (start / m_tick);            // Aux::fill_scalar_blob
            row[4] = (int) ((start + span) / m_tick);
            for (const auto& strip : shape.strips()) {
                if (strip.layer < 2 || strip.layer > 4) continue;
                row[5 + 2 * (strip.layer - 2)] = (int) strip.bounds.first;
                row[6 + 2 * (strip.layer - 2)] = (int) strip.bounds.second;
            }
            // img_eval_all.blobs: yz = (cor * m).sum(1) / max(nc, 1) / 10 -- a sequential sum over the corners
            // as stored in the archive.  Those are the imaging-time corners: ClusterFileSource with
            // restore_corners carries them onto the SimpleBlob.  The reloaded shape re-derives its corners and
            // differs for ~30 % of FD-VD blobs (ClusterArrays loader "fixme"), so it is only the fallback.
            double sy = 0, sz = 0;
            size_t ncorner = 0;
            const auto sblob = std::dynamic_pointer_cast<const Aux::SimpleBlob>(iblob);
            if (m_stored_corners && sblob && !sblob->stored_corners().empty()) {
                for (const auto& c : sblob->stored_corners()) {
                    sy += c.y();
                    sz += c.z();
                }
                ncorner = sblob->stored_corners().size();
                ++m_nstored;
            }
            else {
                for (const auto& c : shape.corners()) {
                    const auto pos = coords.ray_crossing(c.first, c.second);
                    sy += pos[1];
                    sz += pos[2];
                }
                ncorner = shape.corners().size();
            }
            const double nc = std::max<double>((double) ncorner, 1.0);
            const double tmid = (start + 0.5 * span) / 1e3 - m_toff_us;
            row[11] = (m_xw_mm - m_v_mm_per_us * tmid) / 10.0;
            row[12] = sy / nc / 10.0;
            row[13] = sz / nc / 10.0;
            row[14] = (double) iblob->value();
            row[15] = start;
            row[16] = span;
            row[17] = (double) port;
            rows.insert(rows.end(), row, row + NCOL);
        }
    }
    Configuration md;
    md["name"] = "blobs";
    auto tv = std::make_shared<ITensor::vector>();
    tv->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{order, (size_t) NCOL}, rows.data(), md));
    Configuration smd;
    smd["producer"] = "FdvdBlobTable";
    smd["n_stored_corners"] = (Json::UInt64) m_nstored;
    out = std::make_shared<Aux::SimpleTensorSet>(invec[0]->ident(), smd, ITensor::shared_vector(tv));
    log->debug("call={} ident={} blobs={} (stored corners {}) inputs={}", m_count, invec[0]->ident(), order, m_nstored,
               invec.size());
    ++m_count;
    return true;
}
