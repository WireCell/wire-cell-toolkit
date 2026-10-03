#include "WireCellMatch/FdvdLowE.h"

#include "WireCellClus/Facade.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace WireCell;
namespace L = WireCell::Match::FdvdLowE;

// numpy pairwise_sum (numpy/_core/src/umath/loops_utils.h.src, DOUBLE_pairwise_sum, unit stride).  Duplicate of
// flash/src/FdvdOpHitFinder.cxx fdvd_numpy_sum (plugins do not depend on each other).
double L::np_sum(const double* a, size_t n)
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
    return np_sum(a, n2) + np_sum(a + n2, n - n2);
}

// numpy FLOAT_pairwise_sum: the same tree with float accumulators
float L::np_sum_f(const float* a, size_t n)
{
    if (n < 8) {
        float res = 0.0f;
        for (size_t i = 0; i < n; ++i) res += a[i];
        return res;
    }
    if (n <= 128) {
        float r[8];
        for (size_t j = 0; j < 8; ++j) r[j] = a[j];
        size_t i = 8;
        for (; i < n - (n % 8); i += 8) {
            for (size_t j = 0; j < 8; ++j) r[j] += a[i + j];
        }
        float res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    size_t n2 = n / 2;
    n2 -= n2 % 8;
    return np_sum_f(a, n2) + np_sum_f(a + n2, n - n2);
}

std::vector<L::JoinedCluster> L::join_clusters(const Clus::Facade::Grouping& grouping, const double* B, size_t nrow,
                                               size_t ncol, size_t& n_tree_blobs, size_t& n_unjoined)
{
    using Key = std::array<int, 10>;   // anode, face, slice min/max, u/v/w wire min/max
    std::map<Key, int> key2row;
    for (size_t r = 0; r < nrow; ++r) {
        Key k;
        for (int c = 0; c < 10; ++c) k[c] = (int) B[r * ncol + 1 + c];
        key2row[k] = (int) r;
    }
    std::vector<JoinedCluster> out;
    n_tree_blobs = n_unjoined = 0;
    const auto clusters = grouping.children();
    for (size_t ci = 0; ci < clusters.size(); ++ci) {
        JoinedCluster jc;
        jc.index = (int) ci;
        for (const auto* blob : clusters[ci]->children()) {
            const auto wpid = blob->wpid();
            const Key k{wpid.apa(), wpid.face(), blob->slice_index_min(), blob->slice_index_max(),
                        blob->u_wire_index_min(), blob->u_wire_index_max(), blob->v_wire_index_min(),
                        blob->v_wire_index_max(), blob->w_wire_index_min(), blob->w_wire_index_max()};
            ++n_tree_blobs;
            const auto it = key2row.find(k);
            if (it == key2row.end()) {
                ++n_unjoined;
                continue;
            }
            jc.rows.push_back(it->second);
        }
        std::sort(jc.rows.begin(), jc.rows.end(), [&](int a, int b) { return B[a * ncol + kOrder] < B[b * ncol + kOrder]; });
        out.push_back(std::move(jc));
    }
    return out;
}

// numpy/_core/src/multiarray/compiled_base.c arr_interp: binary search, exact-node and last-node short cuts,
// slope * (x - xp[j]) + fp[j]
double L::np_interp(double x, const std::vector<double>& xp, const std::vector<double>& fp, double left, double right)
{
    const size_t n = xp.size();
    if (std::isnan(x)) return x;
    if (x < xp[0]) return left;
    if (x > xp[n - 1]) return right;
    if (x == xp[n - 1]) return fp[n - 1];
    // xp[j] <= x < xp[j+1]
    const size_t j = std::upper_bound(xp.begin(), xp.end(), x) - xp.begin() - 1;
    if (xp[j] == x) return fp[j];
    const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    double res = slope * (x - xp[j]) + fp[j];
    if (std::isnan(res)) {
        res = slope * (x - xp[j + 1]) + fp[j + 1];
        if (std::isnan(res) && fp[j] == fp[j + 1]) res = fp[j];
    }
    return res;
}

// ql_m2m_proto.py shift_table, lines 98-116
bool L::shift_table(const PhotonLibraryModel& lib, const Constants& C, const std::vector<std::array<double, 3>>& P,
                    const std::vector<double>& q, Cluster& c)
{
    const size_t nb = P.size();
    if (nb == 0) return false;
    double pmin = P[0][0], pmax = P[0][0];
    for (const auto& p : P) {
        pmin = std::min(pmin, p[0]);
        pmax = std::max(pmax, p[0]);
    }
    const double t_lo = std::max((C.xc - C.ext_c - pmin) / C.vcm, -C.twin_us);
    const double t_hi = std::min((C.xa + C.ext_a - pmax) / C.vcm, C.twin_us);
    if (t_hi < t_lo) return false;
    const double dt = C.shift_step_cm / C.vcm;
    const size_t n = (size_t) std::ceil((t_hi - t_lo) / dt) + 1;
    c.t_lo = t_lo;
    c.dt = dt;
    c.t_hi = t_hi;
    c.tab.assign(n * NCH, 0.0f);
    std::vector<double> acc(NCH), vis(NCH);
    for (size_t j = 0; j < n; ++j) {
        const double ts = std::min(t_lo + dt * (double) j, t_hi);
        std::fill(acc.begin(), acc.end(), 0.0);
        for (size_t b = 0; b < nb; ++b) {
            const double xt = P[b][0] + C.vcm * ts;
            // numpy's AVX512 exp may differ from libm's in the last bit (doc 16 sec 4.1)
            const double life = std::exp((C.xw_cm - xt) / C.vcm / C.tau_us);
            const double w = q[b] * life;
            lib.visibilities(vis, Point(xt, P[b][1], P[b][2]));
            for (int ch = 0; ch < NCH; ++ch) acc[ch] += vis[ch] * w;   // axis-1 sum: sequential over blobs
        }
        for (int ch = 0; ch < NCH; ++ch) c.tab[j * NCH + ch] = (float) (acc[ch] * C.pe_per_e);
    }
    return true;
}

// ql_m2m_proto.py pred_at, lines 179-187
void L::pred_at(const Cluster& c, double t, double* out)
{
    const size_t n = c.nstep();
    if (n == 1) {
        for (int ch = 0; ch < NCH; ++ch) out[ch] = c.tab[ch];
        return;
    }
    const double u = (t - c.t_lo) / c.dt;
    long k = (long) std::floor(u);
    k = std::clamp<long>(k, 0, std::max<long>((long) n - 2, 0));
    const double f = std::clamp(u - (double) k, 0.0, 1.0);
    const float* a = &c.tab[k * NCH];
    const float* b = &c.tab[(k + 1) * NCH];
    for (int ch = 0; ch < NCH; ++ch) out[ch] = (1 - f) * (double) a[ch] + f * (double) b[ch];
}

// ql_purity.py Resp.__call__ (lines 160-164) and fire (166-167)
double L::Resp::operator()(double p) const
{
    if (p > xc.back()) return p * tail;
    return np_interp(p, xc, y, 0.0, y.back());
}

double L::Resp::fire(double p) const { return np_interp(p, xc, f, 0.0, f.back()); }

L::Calibration L::load_calibration(const std::string& json_path)
{
    const auto path = Persist::resolve(json_path);
    if (path.empty()) raise<ValueError>("FdvdLowE: cannot resolve calibration '%s'", json_path);
    const auto j = Persist::load(path);
    Calibration cal;
    cal.s0 = j["s0"].asDouble();
    cal.rwin_lo = j["rwin"][0].asDouble();
    cal.rwin_hi = j["rwin"][1].asDouble();
    for (const auto& x : j["resp"]["xc"]) cal.resp.xc.push_back(x.asDouble());
    for (const auto& x : j["resp"]["y"]) cal.resp.y.push_back(x.asDouble());
    for (const auto& x : j["resp"]["f"]) cal.resp.f.push_back(x.asDouble());
    cal.resp.tail = j["resp"]["tail"].asDouble();
    for (const auto& x : j["drift"]["mu"]) cal.drift_mu.push_back(x.asDouble());
    for (const auto& x : j["drift"]["dhat"]) cal.drift_dhat.push_back(x.asDouble());
    for (const auto& x : j["drift"]["shat"]) cal.drift_shat.push_back(x.asDouble());
    if (cal.resp.xc.empty() || cal.resp.xc.size() != cal.resp.y.size() || cal.resp.xc.size() != cal.resp.f.size()
        || cal.drift_mu.empty() || cal.drift_mu.size() != cal.drift_dhat.size()
        || cal.drift_mu.size() != cal.drift_shat.size()) {
        raise<ValueError>("FdvdLowE: malformed calibration '%s'", path);
    }
    return cal;
}

// ql_purity.py build_groups (lines 70-88, delta 0) + ql08.py gr_one storage (lines 186-200)
L::Groups L::build_groups(const std::vector<double>& time_ns, const std::vector<double>& pe_rows, const Constants& C)
{
    const size_t nf = time_ns.size();
    std::vector<double> t(nf), pe(nf * NCH), tot(nf);
    for (size_t i = 0; i < nf; ++i) {
        t[i] = time_ns[i] / 1e3;                                             // ql08.reco_arm: of[:, 0] / 1e3
        for (int ch = 0; ch < NCH; ++ch) pe[i * NCH + ch] = (double) (float) pe_rows[i * NCH + ch];   // .astype(float32)
        tot[i] = np_sum(&pe[i * NCH], NCH);
    }
    std::vector<size_t> order_t(nf), by_tot(nf);
    std::iota(order_t.begin(), order_t.end(), 0);
    std::iota(by_tot.begin(), by_tot.end(), 0);
    std::stable_sort(order_t.begin(), order_t.end(), [&](size_t a, size_t b) { return t[a] < t[b]; });
    std::stable_sort(by_tot.begin(), by_tot.end(), [&](size_t a, size_t b) { return -tot[a] < -tot[b]; });
    std::vector<double> ts(nf);
    for (size_t i = 0; i < nf; ++i) ts[i] = t[order_t[i]];
    std::vector<bool> used(nf, false);
    Groups G;
    std::vector<double> gpe(NCH);
    for (size_t f : by_tot) {
        if (used[f]) continue;
        const size_t lo = std::lower_bound(ts.begin(), ts.end(), t[f] - 0.0) - ts.begin();
        const size_t hi = std::upper_bound(ts.begin(), ts.end(), t[f] + 0.0) - ts.begin();
        std::vector<size_t> mem{f};
        if (C.merge_equal_time) {
            for (size_t i = lo; i < hi; ++i) {
                const size_t m = order_t[i];
                if (!used[m] && m != f) mem.push_back(m);
            }
        }
        for (size_t m : mem) used[m] = true;
        // pe[g].sum(0): sequential over members, anchor first
        for (int ch = 0; ch < NCH; ++ch) gpe[ch] = pe[mem[0] * NCH + ch];
        for (size_t k = 1; k < mem.size(); ++k) {
            for (int ch = 0; ch < NCH; ++ch) gpe[ch] += pe[mem[k] * NCH + ch];
        }
        const double gt = t[mem[0]];
        int gnpd = 0;
        for (int ch = 0; ch < NCH; ++ch) gnpd += gpe[ch] >= 1.0;
        if (!(gnpd >= C.pmin_store && std::abs(gt) < C.twin_us + 10)) continue;
        G.t.push_back(gt);
        for (int ch = 0; ch < NCH; ++ch) G.pe.push_back((float) gpe[ch]);
        G.nflash.push_back((int) mem.size());
    }
    return G;
}

namespace {
    // ql_purity.py centroid (line 137): (p @ PD_POS) / max(p.sum(-1), 1e-9); the matmul is BLAS in numpy
    std::array<double, 3> centroid(const double* p, const std::vector<std::array<double, 3>>& pos)
    {
        std::array<double, 3> c{0, 0, 0};
        for (int ch = 0; ch < L::NCH; ++ch) {
            for (int a = 0; a < 3; ++a) c[a] += p[ch] * pos[ch][a];
        }
        const double s = std::max(L::np_sum(p, L::NCH), 1e-9);
        for (int a = 0; a < 3; ++a) c[a] /= s;
        return c;
    }

    // ql_m2m_proto.py ks_dis (lines 196-200)
    double ks_dis(const double* meas, const double* pred)
    {
        const double sa = std::max(L::np_sum(meas, L::NCH), 1e-12);
        const double sb = std::max(L::np_sum(pred, L::NCH), 1e-12);
        double ca = 0, cb = 0, mx = 0;
        for (int ch = 0; ch < L::NCH; ++ch) {
            ca += meas[ch] / sa;
            cb += pred[ch] / sb;
            const double d = std::abs(ca - cb);
            if (ch == 0 || d > mx) mx = d;
        }
        return mx;
    }
}  // namespace

// ql_purity.py features (lines 171-206), tshift = 0
L::Features L::features(const std::vector<Cluster>& cl, const Groups& G, const Calibration& cal, const Constants& C,
                        const std::vector<std::array<double, 3>>& pd_pos)
{
    Features Fe;
    const size_t ng_all = G.t.size();
    std::vector<double> gpe;
    for (size_t i = 0; i < ng_all; ++i) {
        const double gt = G.t[i] + 0.0;
        if (!(std::abs(gt) < C.twin_us)) continue;
        Fe.gidx.push_back((int) i);
        Fe.gt.push_back(gt);
        for (int ch = 0; ch < NCH; ++ch) gpe.push_back((double) G.pe[i * NCH + ch]);
    }
    const size_t ng = Fe.gt.size();
    std::vector<std::array<double, 3>> gcen(ng);
    for (size_t i = 0; i < ng; ++i) {
        const double* p = &gpe[i * NCH];
        Fe.gtot.push_back(np_sum(p, NCH));
        int npd = 0;
        for (int ch = 0; ch < NCH; ++ch) npd += p[ch] >= 1;
        Fe.gnpd.push_back(npd);
        gcen[i] = centroid(p, pd_pos);
    }
    std::vector<double> raw(NCH), pr(NCH), fire(NCH);
    for (size_t k = 0; k < cl.size(); ++k) {
        const auto& c = cl[k];
        for (size_t i = 0; i < ng; ++i) {
            if (!(Fe.gt[i] >= c.t_lo && Fe.gt[i] <= c.t_hi)) continue;
            pred_at(c, Fe.gt[i], raw.data());
            for (int ch = 0; ch < NCH; ++ch) {
                pr[ch] = cal.resp(raw[ch]) * cal.s0;
                fire[ch] = cal.resp.fire(raw[ch]);
            }
            const double npp = np_sum(fire.data(), NCH);
            const double pt = np_sum(pr.data(), NCH);
            const double* meas = &gpe[i * NCH];
            const auto pc = centroid(pr.data(), pd_pos);
            double d2 = 0;
            for (int a = 0; a < 3; ++a) {
                const double d = gcen[i][a] - pc[a];
                d2 += d * d;
            }
            Fe.k.push_back((int) k);
            Fe.g.push_back((int) i);
            Fe.r.push_back((float) (Fe.gtot[i] / std::max(pt, 1e-9)));
            Fe.ks.push_back((float) ks_dis(meas, pr.data()));
            Fe.dc.push_back((float) std::sqrt(d2));
            Fe.npp.push_back((float) npp);
        }
    }
    return Fe;
}

// ql08.py ARMS (lines 69-70)
std::vector<L::Arm> L::default_arms()
{
    return {Arm{"P5E100", 100e3, 5, 100, 0.3, 200.0}, Arm{"P8E0", 100e3, 8, 0, 0.3, 200.0}};
}

// ql10_drift.py drift_mask (lines 77-91)
std::vector<bool> L::drift_mask(const Features& Fe, double veto_k, const std::vector<double>& mu,
                                const std::vector<Cluster>& cl, const Calibration& cal, const Constants& C)
{
    std::vector<bool> ok(Fe.k.size(), true);
    if (veto_k <= 0) return ok;
    for (size_t p = 0; p < ok.size(); ++p) {
        const double mk = mu[Fe.k[p]];
        if (!std::isfinite(mk)) continue;
        const double df = C.x_resp_cm - (cl[Fe.k[p]].x_app + C.vcm * Fe.gt[Fe.g[p]]);
        const double dhat = np_interp(mk, cal.drift_mu, cal.drift_dhat, cal.drift_dhat.front(), cal.drift_dhat.back());
        const double shat = np_interp(mk, cal.drift_mu, cal.drift_shat, cal.drift_shat.front(), cal.drift_shat.back());
        ok[p] = std::abs(df - dhat) <= veto_k * shat;
    }
    return ok;
}

// ql_purity.py decide (lines 209-238), uniq "strict".  NumPy 2 compares the float32 features with the Python
// float cuts in float32 (NEP 50 weak scalars), so the cuts are rounded to float here too.
std::map<int, int> L::decide(const Features& Fe, const std::vector<bool>& mask, const std::vector<Cluster>& cl,
                             const Arm& arm, const Calibration& cal)
{
    float rlo = (float) cal.rwin_lo, rhi = (float) cal.rwin_hi;
    if (arm.rwin_shrink != 1.0) {   // score21.py dec rw 'tight': (rwin[0] * 1.25, rwin[1] / 1.25)
        rlo = (float) (cal.rwin_lo * arm.rwin_shrink);
        rhi = (float) (cal.rwin_hi / arm.rwin_shrink);
    }
    const float ksc = (float) arm.ks_c, dcc = (float) arm.dc_c;
    const float nmin = (float) std::ceil(arm.P / 2.0);
    std::map<int, int> ncl, ngr, best;
    std::vector<std::pair<int, int>> okp;
    for (size_t p = 0; p < Fe.k.size(); ++p) {
        if (!mask[p]) continue;
        const int k = Fe.k[p], g = Fe.g[p];
        const bool ok = cl[k].Q >= arm.qc && Fe.gnpd[g] >= arm.P && Fe.gtot[g] >= arm.E && Fe.r[p] >= rlo
                        && Fe.r[p] <= rhi && Fe.ks[p] <= ksc && Fe.dc[p] <= dcc && Fe.npp[p] >= nmin;
        if (!ok) continue;
        ++ncl[k];
        ++ngr[g];
        okp.emplace_back(k, g);
    }
    std::map<int, int> out;
    for (const auto& [k, g] : okp) {
        if (ncl[k] == 1 && ngr[g] == 1) out[k] = g;
    }
    return out;
}
