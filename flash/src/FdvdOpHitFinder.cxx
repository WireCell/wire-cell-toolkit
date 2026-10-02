#include "WireCellFlash/FdvdOpHitFinder.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include <algorithm>
#include <cmath>
#include <cstring>

WIRECELL_FACTORY(FdvdOpHitFinder, WireCell::Flash::FdvdOpHitFinder,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;

// numpy pairwise_sum (numpy/_core/src/umath/loops_utils.h.src, DOUBLE_pairwise_sum, unit stride):
// below 8 a running sum from 0; up to the 128 block 8 partial sums combined as
// ((r0+r1)+(r2+r3))+((r4+r5)+(r6+r7)) plus the remainder; above, split at n/2 rounded down to a multiple
// of 8.  Checked equal to numpy 2.1.1 sum() on 20000 random float64 arrays (fdvd_sim doc 16 sec 3).
double Flash::fdvd_numpy_sum(const double* a, size_t n)
{
    if (n < 8) {
        double res = 0.0;
        for (size_t i = 0; i < n; ++i) res += a[i];
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (size_t j = 0; j < 8; ++j) r[j] = a[j];
        size_t i = 8;
        for (; i < n - (n % 8); i += 8) {
            for (size_t j = 0; j < 8; ++j) r[j] += a[i + j];
        }
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    size_t n2 = n / 2;
    n2 -= n2 % 8;
    return fdvd_numpy_sum(a, n2) + fdvd_numpy_sum(a + n2, n - n2);
}

// overlay_light.py hitfind (lines 41-75), line for line.
std::vector<Flash::FdvdOpHit> Flash::fdvd_find_ophits(const std::vector<int>& ch, const std::vector<double>& t,
                                                      const std::vector<int64_t>& off,
                                                      const std::vector<double>& adc, const FdvdOpHitParams& par)
{
    const size_t nsn = ch.size();
    const size_t ns = adc.size();
    std::vector<FdvdOpHit> hits;
    if (nsn == 0 || ns == 0) return hits;
    // pedestal: mean of the first nped samples of each snippet, index clamped to the last sample (python :46-48)
    std::vector<double> v(ns);
    std::vector<int64_t> sn(ns);
    for (size_t k = 0; k < nsn; ++k) {
        double sum = 0.0;
        for (int j = 0; j < par.nped; ++j) {
            sum += adc[std::min<int64_t>(off[k] + j, (int64_t) ns - 1)];
        }
        const double ped = sum / par.nped;
        for (int64_t s = off[k]; s < off[k + 1]; ++s) {
            v[s] = adc[s] - ped;
            sn[s] = k;
        }
    }
    // pulses: maximal runs of samples >= threshold2 inside one snippet (python :50-57)
    auto above = [&](size_t s) { return v[s] >= par.threshold2; };
    for (size_t k = 0; k < nsn; ++k) {
        int64_t s = off[k];
        const int64_t s_end = off[k + 1];
        while (s < s_end) {
            if (!above(s)) {
                ++s;
                continue;
            }
            int64_t e = s;
            while (e + 1 < s_end && above(e + 1)) ++e;
            // width = end - start, keep >= min_width (python :58-59)
            if (e - s >= par.min_width) {
                const double* seg = v.data() + s;
                const size_t n = e - s + 1;
                size_t imax = 0;
                for (size_t j = 1; j < n; ++j) {
                    if (seg[j] > seg[imax]) imax = j;   // np.argmax: first maximum
                }
                const double pk = seg[imax];
                if (!(pk < par.threshold || pk < par.hit_threshold)) {   // AlgoSiPM record_hit; HitThreshold
                    const int64_t tmax = s + (int64_t) imax - off[k];
                    const double ts = t[k];
                    const double area = fdvd_numpy_sum(seg, n);
                    FdvdOpHit h;
                    h.channel = ch[k];
                    h.peak_time_us = ts + par.tick_us * (double) tmax;
                    h.start_time_us = ts + par.tick_us * (double) (s - off[k]);
                    h.width_us = par.tick_us * (double) (e - s);
                    h.area = area;
                    h.amplitude = pk;
                    h.pe = area / par.spe_area + par.spe_shift;
                    hits.push_back(h);
                }
            }
            s = e + 1;
        }
    }
    return hits;
}

Flash::FdvdOpHitFinder::FdvdOpHitFinder()
  : Aux::Logger("FdvdOpHitFinder", "flash")
{
}

Flash::FdvdOpHitFinder::~FdvdOpHitFinder() {}

WireCell::Configuration Flash::FdvdOpHitFinder::default_configuration() const
{
    Configuration cfg;
    cfg["nped"] = m_par.nped;
    cfg["threshold"] = m_par.threshold;
    cfg["threshold2"] = m_par.threshold2;
    cfg["min_width"] = m_par.min_width;
    cfg["hit_threshold"] = m_par.hit_threshold;
    cfg["spe_area"] = m_par.spe_area;
    cfg["spe_shift"] = m_par.spe_shift;
    cfg["tick_us"] = m_par.tick_us;
    cfg["max_abs_peak_us"] = m_max_abs_peak_us;
    return cfg;
}

void Flash::FdvdOpHitFinder::configure(const WireCell::Configuration& cfg)
{
    m_par.nped = get(cfg, "nped", m_par.nped);
    m_par.threshold = get(cfg, "threshold", m_par.threshold);
    m_par.threshold2 = get(cfg, "threshold2", m_par.threshold2);
    m_par.min_width = get(cfg, "min_width", m_par.min_width);
    m_par.hit_threshold = get(cfg, "hit_threshold", m_par.hit_threshold);
    m_par.spe_area = get(cfg, "spe_area", m_par.spe_area);
    m_par.spe_shift = get(cfg, "spe_shift", m_par.spe_shift);
    m_par.tick_us = get(cfg, "tick_us", m_par.tick_us);
    m_max_abs_peak_us = get(cfg, "max_abs_peak_us", m_max_abs_peak_us);
    if (m_par.nped < 1 || m_par.spe_area <= 0 || m_par.tick_us <= 0) {
        THROW(ValueError() << errmsg{"FdvdOpHitFinder: nped >= 1, spe_area > 0 and tick_us > 0 required"});
    }
}

namespace {
    ITensor::pointer find_tensor(const ITensorSet::pointer& in, const std::string& name)
    {
        for (const auto& ten : *in->tensors()) {
            if (ten->metadata()["name"].asString() == name) return ten;
        }
        THROW(ValueError() << errmsg{"FdvdOpHitFinder: input tensor set has no tensor named " + name});
    }

    // Any integer or float numpy dtype -> T, element by element (the archive dtypes are the digitiser's).
    template <typename T>
    std::vector<T> as_vector(const ITensor::pointer& ten)
    {
        std::string dt = ten->dtype();
        if (!dt.empty() && (dt[0] == '<' || dt[0] == '|' || dt[0] == '=')) dt = dt.substr(1);
        const size_t n = ten->size() / ten->element_size();
        std::vector<T> out(n);
        const std::byte* p = ten->data();
        auto fill = [&](auto tag) {
            using S = decltype(tag);
            for (size_t i = 0; i < n; ++i) {
                S x;
                std::memcpy(&x, p + i * sizeof(S), sizeof(S));
                out[i] = (T) x;
            }
        };
        if (dt == "i1") fill(int8_t{});
        else if (dt == "u1") fill(uint8_t{});
        else if (dt == "i2") fill(int16_t{});
        else if (dt == "u2") fill(uint16_t{});
        else if (dt == "i4") fill(int32_t{});
        else if (dt == "u4") fill(uint32_t{});
        else if (dt == "i8") fill(int64_t{});
        else if (dt == "f4") fill(float{});
        else if (dt == "f8") fill(double{});
        else {
            THROW(ValueError() << errmsg{"FdvdOpHitFinder: unsupported dtype " + ten->dtype()});
        }
        return out;
    }
}  // namespace

bool Flash::FdvdOpHitFinder::operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    const auto ch = as_vector<int>(find_tensor(in, "ch"));
    const auto t = as_vector<double>(find_tensor(in, "t"));
    const auto off = as_vector<int64_t>(find_tensor(in, "off"));
    const auto adc = as_vector<double>(find_tensor(in, "adc"));
    if (t.size() != ch.size() || off.size() != ch.size() + 1 || (!off.empty() && (size_t) off.back() != adc.size())) {
        THROW(ValueError() << errmsg{"FdvdOpHitFinder: inconsistent snippet tensors"});
    }
    const auto hits = fdvd_find_ophits(ch, t, off, adc, m_par);

    // ophits_to_tensor.py: drop |peak| > max_abs_peak_us, times x 1000 to ns
    std::vector<double> flat;
    flat.reserve(hits.size() * 9);
    size_t nkeep = 0, ndrop = 0;
    for (const auto& h : hits) {
        if (std::abs(h.peak_time_us) > m_max_abs_peak_us) {
            ++ndrop;
            continue;
        }
        flat.push_back((double) h.channel);
        flat.push_back(h.peak_time_us * 1000.0);
        flat.push_back(h.width_us * 1000.0);
        flat.push_back(h.area);
        flat.push_back(h.amplitude);
        flat.push_back(h.pe);
        flat.push_back(h.start_time_us * 1000.0);
        flat.push_back(-1.0);
        flat.push_back(0.0);
        ++nkeep;
    }
    Configuration tmd;
    tmd["name"] = "ophits";
    auto tensors = std::make_shared<ITensor::vector>();
    tensors->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{nkeep, (size_t) 9}, flat.data(), tmd));
    Configuration md = in->metadata();
    md["producer"] = "FdvdOpHitFinder";
    md["n_dropped"] = (Json::UInt64) ndrop;
    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), md, ITensor::shared_vector(tensors));
    log->debug("call={} ident={} snippets={} hits={} dropped={}", m_count, in->ident(), ch.size(), nkeep, ndrop);
    ++m_count;
    return true;
}
