#include "WireCellMatch/FdvdDriftRegressor3View.h"
#include "WireCellMatch/FdvdLowE.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
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

FdvdDriftRegressor3View::FdvdDriftRegressor3View(const WireCell::Configuration& cfg)
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
        raise<ValueError>("FdvdDriftRegressor3View: frames must be a path pattern with %%d for the anode, got '%s'", m_frames);
    }
    m_fwd = Factory::find_tn<ITensorForward>(m_forward);
}

// apply35.py collect_one: np.abs(frame) of the tagged frame, float32, every row
std::vector<float> FdvdDriftRegressor3View::load_rows(const std::string& path, const std::string& frame_tag, size_t& nrow,
                                                      size_t& ntick, int& ch0)
{
    boost::iostreams::filtering_istream in;
    Stream::input_filters(in, path);
    if (in.size() < 1) raise<IOError>("FdvdDriftRegressor3View: cannot open %s", path);
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
            if (!pigenc::eigen::load(pig, frame)) raise<IOError>("FdvdDriftRegressor3View: bad frame %s in %s", fname, path);
            have_frame = true;
        }
        else if (tagged && fname.rfind("channels_", 0) == 0) {
            pigenc::stl::load(pig, channels);
        }
    }
    if (!have_frame || channels.empty() || (size_t) frame.rows() != channels.size()) {
        raise<IOError>("FdvdDriftRegressor3View: no usable '%s' frame in %s", frame_tag, path);
    }
    for (size_t k = 1; k < channels.size(); ++k) {
        if (channels[k] != channels[k - 1] + 1) raise<IOError>("FdvdDriftRegressor3View: non-contiguous channels in %s", path);
    }
    nrow = frame.rows();
    ntick = frame.cols();
    ch0 = channels.front();
    std::vector<float> img(nrow * ntick);
    for (size_t r = 0; r < nrow; ++r) {
        for (size_t t = 0; t < ntick; ++t) img[r * ntick + t] = std::abs(frame(r, t));
    }
    return img;
}

// apply35.py collect_one, one view: footprint mask, masked image, charge centroid, cut (the numpy order of
// FdvdDriftRegressor::compute: float32 pairwise sums, float64 centroid, round half to even)
bool FdvdDriftRegressor3View::view_crop(const float* img, int nrow, size_t ntick, const std::vector<box_t>& boxes, int dw,
                                        int dt, int nch, int ntk, float* out, long& c0, long& t0)
{
    std::fill(out, out + (size_t) nch * ntk, 0.0f);
    c0 = 0;
    t0 = 0;
    std::vector<char> m((size_t) nrow * ntick, 0);
    for (const auto& b : boxes) {
        const long r0 = std::max(b[0] - dw, 0L), r1 = std::min(b[1] + dw, (long) nrow);
        const long a0 = std::max(b[2] - dt, 0L), a1 = std::min(b[3] + dt, (long) ntick);
        for (long rr = r0; rr < r1; ++rr) {
            for (long cc = a0; cc < a1; ++cc) m[rr * ntick + cc] = 1;
        }
    }
    std::vector<float> mi((size_t) nrow * ntick);
    for (size_t p = 0; p < mi.size(); ++p) mi[p] = m[p] ? img[p] : 0.0f;
    const float w = L::np_sum_f(mi.data(), mi.size());
    if (w <= 0) return false;
    std::vector<double> rw(nrow), cwv(ntick);
    for (int rr = 0; rr < nrow; ++rr) rw[rr] = (double) L::np_sum_f(&mi[rr * ntick], ntick) * (double) rr;
    std::vector<float> col(mi.begin(), mi.begin() + ntick);
    for (int rr = 1; rr < nrow; ++rr) {
        for (size_t cc = 0; cc < ntick; ++cc) col[cc] += mi[rr * ntick + cc];
    }
    for (size_t cc = 0; cc < ntick; ++cc) cwv[cc] = (double) col[cc] * (double) cc;
    const double cw = L::np_sum(rw.data(), rw.size()) / (double) w;
    const double ct = L::np_sum(cwv.data(), cwv.size()) / (double) w;
    c0 = (long) std::nearbyint(cw) - nch / 2;
    t0 = (long) std::nearbyint(ct) - ntk / 2;
    for (long rr = std::max(c0, 0L); rr < std::min(c0 + nch, (long) nrow); ++rr) {
        for (long cc = std::max(t0, 0L); cc < std::min(t0 + ntk, (long) ntick); ++cc) {
            out[(rr - c0) * ntk + (cc - t0)] = mi[rr * ntick + cc];
        }
    }
    return true;
}

// score_val_uvw.py fusion_input: x[k][:, lo:hi] = src[:, lo - dt:hi - dt]
void FdvdDriftRegressor3View::shift_ticks(const float* src, long d, int nch, int ntk, float* dst)
{
    std::fill(dst, dst + (size_t) nch * ntk, 0.0f);
    const long lo = std::max(0L, d), hi = std::min((long) ntk, (long) ntk + d);
    if (hi <= lo) return;
    for (int r = 0; r < nch; ++r) {
        std::copy(src + (size_t) r * ntk + (lo - d), src + (size_t) r * ntk + (hi - d), dst + (size_t) r * ntk + lo);
    }
}

FdvdDriftRegressor3View::Result FdvdDriftRegressor3View::compute(const Clus::Facade::Grouping& grouping, const double* B,
                                                                 size_t nbt, size_t nc, int ident) const
{
    size_t ntree = 0, nmiss = 0;
    const auto joined = L::join_clusters(grouping, B, nbt, nc, ntree, nmiss);

    // selection and dominant CRM: FdvdDriftRegressor::compute
    struct Sel {
        int rep;
        int crm;
        std::vector<int> rows;   // blobs in the dominant CRM, python order
    };
    std::map<int, std::vector<Sel>> by_crm;
    for (const auto& jc : joined) {
        if (jc.rows.empty()) continue;
        double qsig = 0;
        std::vector<double> pq;
        std::map<int, double> qcrm;
        for (int r : jc.rows) {
            const double q = B[r * nc + L::kQ];
            qsig += q;
            pq.push_back(std::max(q, 0.0));
            qcrm[(int) B[r * nc + L::kAnode]] += pq.back();
        }
        const double qpos = L::np_sum(pq.data(), pq.size());
        if (!(qsig >= m_qmin_all) || qpos <= 0 || !(qpos >= m_qmin)) continue;
        int crm = -1;
        double best = -1;
        for (const auto& [c, q] : qcrm) {
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
    R.ncrm = by_crm.size();
    const size_t csz = (size_t) m_nch * m_ntk;
    // blob-table columns of the wire range per view: FdvdBlobTable row[5 + 2 * (layer - 2)], row[6 + 2 * (layer - 2)]
    const int wcol[3] = {5, 7, L::kWmin};
    for (const auto& [crm, sels] : by_crm) {
        size_t nrow = 0, ntick = 0;
        int ch0 = 0;
        const auto img = load_rows(String::format(m_frames, crm), m_frame_tag, nrow, ntick, ch0);
        if ((int) nrow < m_row0[2] + m_nrow[2]) raise<IOError>("FdvdDriftRegressor3View: frame of CRM %d has %d rows", crm, (int) nrow);
        for (size_t j0 = 0; j0 < sels.size(); j0 += m_chunk) {
            const size_t j1 = std::min(sels.size(), j0 + (size_t) m_chunk), nj = j1 - j0;
            std::vector<float> own(nj * 3 * csz, 0.0f);     // the three own crops
            std::vector<long> org(nj * 6, 0);               // c0, t0 per view
            std::vector<char> empty(nj * 3, 0);
            for (size_t j = j0; j < j1; ++j) {
                long* o = &org[(j - j0) * 6];
                for (int p = 0; p < 3; ++p) {
                    std::vector<box_t> boxes;
                    for (int r : sels[j].rows) {
                        const double* b = B + r * nc;
                        const long t0 = (long) std::floor(b[L::kStart] / m_tick);
                        boxes.push_back({(long) (int) b[wcol[p]], (long) (int) b[wcol[p] + 1], t0,
                                         t0 + (long) std::ceil(b[L::kSpan] / m_tick)});
                    }
                    const bool ok = view_crop(&img[(size_t) m_row0[p] * ntick], m_nrow[p], ntick, boxes, m_dw, m_dt, m_nch,
                                              m_ntk, &own[((j - j0) * 3 + p) * csz], o[2 * p], o[2 * p + 1]);
                    empty[(j - j0) * 3 + p] = !ok;
                }
                for (int p = 0; p < 2; ++p) {               // an empty induction view takes the W origin
                    if (empty[(j - j0) * 3 + p] && !empty[(j - j0) * 3 + 2]) {
                        o[2 * p] = o[4];
                        o[2 * p + 1] = o[5];
                    }
                }
            }
            for (size_t k0 = j0; k0 < j1; k0 += m_batch) {
                const size_t k1 = std::min(j1, k0 + (size_t) m_batch), nb = k1 - k0;
                std::vector<float> x(nb * 3 * csz), meta(nb * 4);
                std::vector<float> sh(csz);
                for (size_t n = 0; n < nb; ++n) {
                    const size_t jj = k0 - j0 + n;
                    const long* o = &org[jj * 6];
                    for (int p = 0; p < 3; ++p) {
                        const float* src = &own[(jj * 3 + p) * csz];
                        shift_ticks(src, o[2 * p + 1] - o[5], m_nch, m_ntk, sh.data());
                        float* dst = &x[(n * 3 + p) * csz];
                        for (size_t q = 0; q < csz; ++q) dst[q] = std::log1p(sh[q]) / 5.0f;
                        meta[n * 4 + 1 + p] = (float) (ch0 + m_row0[p] + o[2 * p]);
                    }
                    meta[n * 4] = (float) crm;
                }
                Configuration md, mm;
                md["name"] = "crops";
                mm["name"] = "meta";
                auto tv = std::make_shared<ITensor::vector>();
                tv->push_back(std::make_shared<Aux::SimpleTensor>(
                    ITensor::shape_t{nb, (size_t) 3, (size_t) m_nch, (size_t) m_ntk}, x.data(), md));
                tv->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{nb, (size_t) 4}, meta.data(), mm));
                auto res = m_fwd->forward(std::make_shared<Aux::SimpleTensorSet>(ident, Configuration{},
                                                                                 ITensor::shared_vector(tv)));
                if (!res || res->tensors()->size() < 2) raise<ValueError>("FdvdDriftRegressor3View: model returned no (mu, logvar)");
                const float* mu = (const float*) res->tensors()->at(0)->data();
                const float* lv = (const float*) res->tensors()->at(1)->data();
                for (size_t n = 0; n < nb; ++n) {
                    const float mcm = mu[n] * 100.0f;
                    const float sig = std::exp(0.5f * std::clamp(lv[n], -7.0f, 7.0f)) * 100.0f;
                    R.rows.insert(R.rows.end(), {(double) sels[k0 + n].rep, (double) mcm, (double) sig});
                }
            }
            if (m_dump_crops) {
                R.crops.insert(R.crops.end(), own.begin(), own.end());
                for (size_t jj = 0; jj < nj; ++jj) {
                    R.origins.push_back((double) crm);
                    for (int p = 0; p < 3; ++p) {
                        R.origins.push_back((double) org[jj * 6 + 2 * p]);
                        R.origins.push_back(p < 2 && empty[jj * 3 + p] ? -1e9 : (double) org[jj * 6 + 2 * p + 1]);
                    }
                }
            }
            R.ncrop += nj;
        }
    }

    SPDLOG_LOGGER_DEBUG(m_log, "FdvdDriftRegressor3View ident={} crops {} in {} CRMs", ident, R.ncrop, by_crm.size());
    return R;
}
