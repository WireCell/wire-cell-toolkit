#include "WireCellMatch/FdvdDriftRegressor.h"
#include "WireCellMatch/FdvdLowE.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellClus/Facade.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/custard/pigenc_eigen.hpp"
#include "WireCellUtil/custard/pigenc_stl.hpp"

#include <boost/iostreams/filtering_stream.hpp>

#include <algorithm>
#include <cmath>
#include <map>

using namespace WireCell;
using namespace WireCell::Match;
namespace L = WireCell::Match::FdvdLowE;

FdvdDriftRegressor::FdvdDriftRegressor(const WireCell::Configuration& cfg)
  : m_log(Log::logger("match"))
{
    m_forward = get(cfg, "forward", m_forward);
    m_frames = get(cfg, "frames", m_frames);
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);
    m_qmin_all = get(cfg, "qmin_all", m_qmin_all);
    m_qmin = get(cfg, "qmin", m_qmin);
    m_tick = get(cfg, "tick", m_tick);
    m_chunk = get(cfg, "chunk", m_chunk);
    m_batch = get(cfg, "batch", m_batch);
    m_dump_crops = get(cfg, "dump_crops", m_dump_crops);
    if (m_frames.find('%') == std::string::npos) {
        raise<ValueError>("FdvdDriftRegressor: frames must be a path pattern with %%d for the anode, got '%s'", m_frames);
    }
    m_fwd = Factory::find_tn<ITensorForward>(m_forward);
}

// dl_drift.py load_gauss (lines 175-183): |frame| rows w0 .. w0+nw of the tagged frame, float32
std::vector<float> FdvdDriftRegressor::load_w_rows(const std::string& path, const std::string& frame_tag, int w0, int nw,
                                                   size_t& ntick)
{
    boost::iostreams::filtering_istream in;
    Stream::input_filters(in, path);
    if (in.size() < 1) raise<IOError>("FdvdDriftRegressor: cannot open %s", path);
    using array_t = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    array_t frame;
    std::vector<int> channels;
    bool have_frame = false;
    while (true) {
        std::string fname;
        size_t fsize = 0;
        custard::read(in, fname, fsize);
        if (fsize == 0 || !in) break;
        pigenc::File pig;
        pig.read(in);
        const bool tagged = fname.find("_" + frame_tag) != std::string::npos;
        if (tagged && fname.rfind("frame_", 0) == 0) {
            if (!pigenc::eigen::load(pig, frame)) raise<IOError>("FdvdDriftRegressor: bad frame %s in %s", fname, path);
            have_frame = true;
        }
        else if (tagged && fname.rfind("channels_", 0) == 0) {
            pigenc::stl::load(pig, channels);
        }
    }
    if (!have_frame || frame.rows() < w0 + nw) raise<IOError>("FdvdDriftRegressor: no usable '%s' frame in %s", frame_tag, path);
    for (size_t k = 1; k < channels.size(); ++k) {
        if (channels[k] != channels[k - 1] + 1) raise<IOError>("FdvdDriftRegressor: non-contiguous channels in %s", path);
    }
    ntick = frame.cols();
    std::vector<float> img((size_t) nw * ntick);
    for (int r = 0; r < nw; ++r) {
        for (size_t t = 0; t < ntick; ++t) img[r * ntick + t] = std::abs(frame(w0 + r, t));
    }
    return img;
}

FdvdDriftRegressor::Result FdvdDriftRegressor::compute(const Clus::Facade::Grouping& grouping, const double* B,
                                                       size_t nbt, size_t nc, int ident) const
{
    size_t ntree = 0, nmiss = 0;
    const auto joined = L::join_clusters(grouping, B, nbt, nc, ntree, nmiss);

    // dl_drift.collect_one (lines 120-148) + ovl13 Chain.drift (line 371): selection, dominant CRM, footprint
    struct Sel {
        int rep;
        int crm;
        std::vector<int> rows;   // blobs in the dominant CRM, python order
    };
    std::map<int, std::vector<Sel>> by_crm;   // sorted by CRM, as crops_event
    for (const auto& jc : joined) {
        if (jc.rows.empty()) continue;
        double qsig = 0;
        std::vector<double> pq;
        std::map<int, double> qcrm;
        for (int r : jc.rows) {
            const double q = B[r * nc + L::kQ];
            qsig += q;
            pq.push_back(std::max(q, 0.0));
            qcrm[(int) B[r * nc + L::kAnode]] += pq.back();   // np.bincount: sequential
        }
        const double qpos = L::np_sum(pq.data(), pq.size());
        if (!(qsig >= m_qmin_all) || qpos <= 0 || !(qpos >= m_qmin)) continue;
        int crm = -1;
        double best = -1;
        for (const auto& [c, q] : qcrm) {   // argmax: first (lowest CRM) maximum
            if (q > best) {
                best = q;
                crm = c;
            }
        }
        Sel s{(int) B[jc.rows.front() * nc + L::kOrder], crm, {}};
        for (int r : jc.rows) {
            if ((int) B[r * nc + L::kAnode] == crm) s.rows.push_back(r);
        }
        by_crm[crm].push_back(std::move(s));
    }

    Result R;
    auto& drift = R.rows;
    auto& crops = R.crops;
    size_t& ncrop = R.ncrop;
    R.ncrm = by_crm.size();
    const size_t csz = (size_t) m_nch * m_ntk;
    for (const auto& [crm, sels] : by_crm) {
        size_t ntick = 0;
        const auto img = load_w_rows(String::format(m_frames, crm), m_frame_tag, m_w0, m_nw, ntick);
        for (size_t j0 = 0; j0 < sels.size(); j0 += m_chunk) {
            const size_t j1 = std::min(sels.size(), j0 + (size_t) m_chunk);
            std::vector<float> own((j1 - j0) * csz, 0.0f);
            std::vector<float> mi((size_t) m_nw * ntick);
            for (size_t j = j0; j < j1; ++j) {
                // footprint_mask (lines 186-191) and the masked image np.where(m, img, 0.0) (float32)
                std::vector<char> m((size_t) m_nw * ntick, 0);
                for (int r : sels[j].rows) {
                    const double* b = B + r * nc;
                    const int w0 = (int) b[L::kWmin], w1 = (int) b[L::kWmax];
                    const long t0 = (long) std::floor(b[L::kStart] / m_tick);
                    const long t1 = t0 + (long) std::ceil(b[L::kSpan] / m_tick);
                    const int r0 = std::max(w0 - m_dw, 0), r1 = std::min(w1 + m_dw, m_nw);
                    const long c0 = std::max(t0 - m_dt, 0L), c1 = std::min(t1 + m_dt, (long) ntick);
                    for (int rr = r0; rr < r1; ++rr) {
                        for (long cc = c0; cc < c1; ++cc) m[rr * ntick + cc] = 1;
                    }
                }
                for (size_t p = 0; p < mi.size(); ++p) mi[p] = m[p] ? img[p] : 0.0f;
                const float w = L::np_sum_f(mi.data(), mi.size());
                if (w <= 0) continue;   // crop stays zero, still inferred (python)
                // charge centroid: (mi.sum(1) * arange).sum() / w and (mi.sum(0) * arange).sum() / w
                std::vector<double> rw(m_nw), cwv(ntick);
                for (int rr = 0; rr < m_nw; ++rr) rw[rr] = (double) L::np_sum_f(&mi[rr * ntick], ntick) * (double) rr;
                std::vector<float> col(mi.begin(), mi.begin() + ntick);
                for (int rr = 1; rr < m_nw; ++rr) {
                    for (size_t cc = 0; cc < ntick; ++cc) col[cc] += mi[rr * ntick + cc];
                }
                for (size_t cc = 0; cc < ntick; ++cc) cwv[cc] = (double) col[cc] * (double) cc;
                const double cw = L::np_sum(rw.data(), rw.size()) / (double) w;
                const double ct = L::np_sum(cwv.data(), cwv.size()) / (double) w;
                // cut (lines 194-203): int(round(.)) is round half to even
                const long c0 = (long) std::nearbyint(cw) - m_nch / 2, t0 = (long) std::nearbyint(ct) - m_ntk / 2;
                float* o = &own[(j - j0) * csz];
                for (long rr = std::max(c0, 0L); rr < std::min(c0 + m_nch, (long) m_nw); ++rr) {
                    for (long cc = std::max(t0, 0L); cc < std::min(t0 + m_ntk, (long) ntick); ++cc) {
                        o[(rr - c0) * m_ntk + (cc - t0)] = mi[rr * ntick + cc];
                    }
                }
            }
            // dl_drift.predict (lines 233-246) in batches
            for (size_t k0 = j0; k0 < j1; k0 += m_batch) {
                const size_t k1 = std::min(j1, k0 + (size_t) m_batch), nb = k1 - k0;
                std::vector<float> x(nb * csz);
                for (size_t p = 0; p < x.size(); ++p) x[p] = std::log1p(own[(k0 - j0) * csz + p]) / 5.0f;
                Configuration md;
                md["name"] = "crops";
                auto tv = std::make_shared<ITensor::vector>();
                tv->push_back(std::make_shared<Aux::SimpleTensor>(
                    ITensor::shape_t{nb, (size_t) 1, (size_t) m_nch, (size_t) m_ntk}, x.data(), md));
                auto res = m_fwd->forward(std::make_shared<Aux::SimpleTensorSet>(ident, Configuration{},
                                                                                 ITensor::shared_vector(tv)));
                if (!res || res->tensors()->size() < 2) raise<ValueError>("FdvdDriftRegressor: model returned no (mu, logvar)");
                const float* mu = (const float*) res->tensors()->at(0)->data();
                const float* lv = (const float*) res->tensors()->at(1)->data();
                for (size_t n = 0; n < nb; ++n) {
                    const float mcm = mu[n] * 100.0f;
                    const float sig = std::exp(0.5f * std::clamp(lv[n], -7.0f, 7.0f)) * 100.0f;
                    drift.insert(drift.end(), {(double) sels[k0 + n].rep, (double) mcm, (double) sig});
                }
            }
            if (m_dump_crops) crops.insert(crops.end(), own.begin(), own.end());
            ncrop += j1 - j0;
        }
    }

    SPDLOG_LOGGER_DEBUG(m_log, "FdvdDriftRegressor ident={} crops {} in {} CRMs", ident, ncrop, by_crm.size());
    return R;
}
