#include "WireCellFlash/FdvdAdjOpFlashFinder.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>

WIRECELL_FACTORY(FdvdAdjOpFlashFinder, WireCell::Flash::FdvdAdjOpFlashFinder,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;

// AdjOpHitsUtils.cc GetOpHitPlane (:559-589) with CheckPlane (:602-613): open interval of +-buffer
// around each X-ARAPUCA plane, tested in this order (adjophits.py opdet_planes).
int Flash::fdvd_adj_plane(double x, double y, double z, const FdvdAdjParams& par)
{
    const double b = par.plane_buffer;
    if (par.xa_cathode_x - b < x && x < par.xa_cathode_x + b) return 0;
    if (par.xa_membrane_y - b < y && y < par.xa_membrane_y + b) return 1;
    if (-par.xa_membrane_y - b < y && y < -par.xa_membrane_y + b) return 2;
    if (par.xa_final_cap_z - b < z && z < par.xa_final_cap_z + b) return 3;
    if (par.xa_start_cap_z - b < z && z < par.xa_start_cap_z + b) return 4;
    return -1;
}

// AdjOpHitsUtils.cc CalcAdjOpHits (:262-494), the control flow of adjophits.py calc_adj_ophits / adj21.py calc.
std::vector<Flash::FdvdAdjCluster> Flash::fdvd_calc_adj_ophits(const std::vector<double>& t_ns,
                                                               const std::vector<double>& pe,
                                                               const std::vector<int>& opdet,
                                                               const std::vector<std::vector<double>>& dist,
                                                               const std::vector<int>& plane, const FdvdAdjParams& par)
{
    const size_t n = t_ns.size();
    std::vector<size_t> sidx(n);
    std::iota(sidx.begin(), sidx.end(), 0);
    std::stable_sort(sidx.begin(), sidx.end(), [&](size_t a, size_t b) { return t_ns[a] < t_ns[b]; });   // :296 stable_sort
    const double wmax = par.max_time_us * 1000.0;
    const double wmin = par.min_time_us * 1000.0;
    std::vector<char> claimed(n, 0);
    std::vector<FdvdAdjCluster> clusters;
    std::vector<size_t> members;
    for (size_t pos = 0; pos < n; ++pos) {
        const size_t i = sidx[pos];
        if (claimed[i] || pe[i] < par.trigger_pe) continue;          // :301-312 seed
        bool main = true;
        claimed[i] = 1;
        members.assign(1, i);
        const auto& Di = dist[opdet[i]];
        const double ti = t_ns[i];
        const int pi = plane[opdet[i]];
        // forward scan (:319-395)
        for (size_t q = pos + 1; q < n; ++q) {
            const size_t j = sidx[q];
            if (std::abs(t_ns[j] - ti) > wmax) break;
            if (pe[j] < par.pe || pi != plane[opdet[j]] || (claimed[j] && !par.hit_duplicates)) continue;
            if (Di[opdet[j]] < par.radius_cm) {
                if (pe[j] > pe[i]) {                                    // a brighter neighbour demotes the seed
                    main = false;
                    for (size_t m : members) claimed[m] = 0;
                    break;
                }
                members.push_back(j);
                claimed[j] = 1;
            }
        }
        if (!main) {
            claimed[i] = 0;
            continue;
        }
        // backward scan (:403-479): `it4 != begin` never reaches sorted position 0
        if (pos != 0) {
            for (size_t q = pos - 1; q > 0; --q) {
                const size_t j = sidx[q];
                if (std::abs(ti - t_ns[j]) > wmin) break;
                if (pe[j] < par.pe || pi != plane[opdet[j]] || (claimed[j] && !par.hit_duplicates)) continue;
                if (Di[opdet[j]] < par.radius_cm) {
                    if (pe[j] > pe[i]) {
                        main = false;
                        for (size_t m : members) claimed[m] = 0;
                        break;
                    }
                    members.push_back(j);
                    claimed[j] = 1;
                }
            }
        }
        if (main && (int) members.size() >= par.nhit) {
            clusters.push_back(members);
        }
        else if (!main) {
            claimed[i] = 0;
        }
    }
    return clusters;
}

Flash::FdvdAdjOpFlashFinder::FdvdAdjOpFlashFinder()
  : Aux::Logger("FdvdAdjOpFlashFinder", "flash")
{
}

Flash::FdvdAdjOpFlashFinder::~FdvdAdjOpFlashFinder() {}

WireCell::Configuration Flash::FdvdAdjOpFlashFinder::default_configuration() const
{
    Configuration cfg;
    cfg["nchan"] = m_nchan;
    cfg["geom_file"] = m_geom_file;
    cfg["channel_map_file"] = m_channel_map_file;
    cfg["time_var"] = m_par.start_time ? "StartTime" : "PeakTime";
    cfg["min_time_us"] = m_par.min_time_us;
    cfg["max_time_us"] = m_par.max_time_us;
    cfg["radius_cm"] = m_par.radius_cm;
    cfg["nhit"] = m_par.nhit;
    cfg["pe"] = m_par.pe;
    cfg["trigger_pe"] = m_par.trigger_pe;
    cfg["hot_threshold"] = m_par.hot_threshold;
    cfg["hit_duplicates"] = m_par.hit_duplicates;
    cfg["xa_cathode_x"] = m_par.xa_cathode_x;
    cfg["xa_membrane_y"] = m_par.xa_membrane_y;
    cfg["xa_final_cap_z"] = m_par.xa_final_cap_z;
    cfg["xa_start_cap_z"] = m_par.xa_start_cap_z;
    cfg["plane_buffer"] = m_par.plane_buffer;
    return cfg;
}

void Flash::FdvdAdjOpFlashFinder::configure(const WireCell::Configuration& cfg)
{
    m_nchan = get(cfg, "nchan", m_nchan);
    m_geom_file = get(cfg, "geom_file", m_geom_file);
    m_channel_map_file = get(cfg, "channel_map_file", m_channel_map_file);
    const std::string tv = get<std::string>(cfg, "time_var", m_par.start_time ? "StartTime" : "PeakTime");
    if (tv != "StartTime" && tv != "PeakTime") {
        THROW(ValueError() << errmsg{"FdvdAdjOpFlashFinder: time_var must be StartTime or PeakTime"});
    }
    m_par.start_time = tv == "StartTime";
    m_par.min_time_us = get(cfg, "min_time_us", m_par.min_time_us);
    m_par.max_time_us = get(cfg, "max_time_us", m_par.max_time_us);
    m_par.radius_cm = get(cfg, "radius_cm", m_par.radius_cm);
    m_par.nhit = get(cfg, "nhit", m_par.nhit);
    m_par.pe = get(cfg, "pe", m_par.pe);
    m_par.trigger_pe = get(cfg, "trigger_pe", m_par.trigger_pe);
    m_par.hot_threshold = get(cfg, "hot_threshold", m_par.hot_threshold);
    m_par.hit_duplicates = get(cfg, "hit_duplicates", m_par.hit_duplicates);
    m_par.xa_cathode_x = get(cfg, "xa_cathode_x", m_par.xa_cathode_x);
    m_par.xa_membrane_y = get(cfg, "xa_membrane_y", m_par.xa_membrane_y);
    m_par.xa_final_cap_z = get(cfg, "xa_final_cap_z", m_par.xa_final_cap_z);
    m_par.xa_start_cap_z = get(cfg, "xa_start_cap_z", m_par.xa_start_cap_z);
    m_par.plane_buffer = get(cfg, "plane_buffer", m_par.plane_buffer);
    if (m_geom_file.empty() || m_channel_map_file.empty()) {
        THROW(ValueError() << errmsg{"FdvdAdjOpFlashFinder: geom_file and channel_map_file are required"});
    }
    m_chmap.clear();
    const auto jmap = Persist::load(m_channel_map_file);      // keep the loaded value alive over the loop
    for (const auto& jc : jmap["channels"]) {
        m_chmap[jc["opch"].asInt()] = jc["opdet"].asInt();
    }
    // OpDet centres: the json is in mm; adjophits.py opdet_centres() divides by 10
    m_pos.assign(m_nchan, {0.0, 0.0, 0.0});
    std::vector<char> seen(m_nchan, 0);
    const auto jgeom = Persist::load(m_geom_file);
    for (const auto& jod : jgeom["opdets"]) {
        const int od = jod["opdet"].asInt();
        if (od < 0 || od >= m_nchan) continue;
        m_pos[od] = {jod["x"].asDouble() / 10.0, jod["y"].asDouble() / 10.0, jod["z"].asDouble() / 10.0};
        seen[od] = 1;
    }
    if (std::count(seen.begin(), seen.end(), 1) != m_nchan) {
        THROW(ValueError() << errmsg{"FdvdAdjOpFlashFinder: geom_file " + m_geom_file + " holds " +
                                     std::to_string(std::count(seen.begin(), seen.end(), 1)) + " of " +
                                     std::to_string(m_nchan) + " OpDets"});
    }
    m_dist.assign(m_nchan, std::vector<double>(m_nchan, 0.0));
    m_plane.assign(m_nchan, -1);
    for (int a = 0; a < m_nchan; ++a) {
        m_plane[a] = fdvd_adj_plane(m_pos[a][0], m_pos[a][1], m_pos[a][2], m_par);
        for (int b = 0; b < m_nchan; ++b) {
            const double dx = m_pos[a][0] - m_pos[b][0], dy = m_pos[a][1] - m_pos[b][1], dz = m_pos[a][2] - m_pos[b][2];
            m_dist[a][b] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
}

namespace {
    ITensor::pointer find_tensor(const ITensorSet::pointer& in, const std::string& name)
    {
        for (const auto& ten : *in->tensors()) {
            if (ten->metadata()["name"].asString() == name) return ten;
        }
        THROW(ValueError() << errmsg{"FdvdAdjOpFlashFinder: input tensor set has no tensor named " + name});
    }

    ITensor::pointer make_tensor(const std::string& name, const std::vector<double>& data, size_t ncol)
    {
        Configuration md;
        md["name"] = name;
        const size_t nrow = ncol ? data.size() / ncol : 0;
        return std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{nrow, ncol}, data.data(), md);
    }
}  // namespace

bool Flash::FdvdAdjOpFlashFinder::operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    const auto ten = find_tensor(in, "ophits");
    const auto shape = ten->shape();
    if (shape.size() != 2 || shape[1] < 7 || ten->element_size() != sizeof(double)) {
        THROW(ValueError() << errmsg{"FdvdAdjOpFlashFinder: ophits must be f8 [nhit, >= 7]"});
    }
    const size_t nh = shape[0], nc = shape[1];
    std::vector<double> H(nh * nc);
    if (nh) std::memcpy(H.data(), ten->data(), nh * nc * sizeof(double));
    std::vector<double> t(nh), pe(nh);
    std::vector<int> od(nh);
    for (size_t i = 0; i < nh; ++i) {
        const double* r = &H[i * nc];
        t[i] = m_par.start_time ? r[6] : r[1];
        pe[i] = r[5];
        const auto it = m_chmap.find((int) r[0]);
        if (it == m_chmap.end() || it->second < 0 || it->second >= m_nchan) {
            THROW(ValueError() << errmsg{"FdvdAdjOpFlashFinder: OpChannel without an OpDet: " + std::to_string((int) r[0])});
        }
        od[i] = it->second;
    }
    const auto clusters = fdvd_calc_adj_ophits(t, pe, od, m_dist, m_plane, m_par);

    const size_t nf = clusters.size(), mcol = 1 + (size_t) m_nchan;
    std::vector<double> flash(nf * mcol, 0.0), summary(nf * 8, 0.0), member;
    for (size_t f = 0; f < nf; ++f) {
        const auto& u = clusters[f];
        size_t best = 0;                                              // MakeFlashVector: max-PE member, strict >
        for (size_t q = 1; q < u.size(); ++q) {
            if (pe[u[q]] > pe[u[best]]) best = q;
        }
        double* row = &flash[f * mcol];
        row[0] = t[u[best]];
        double tot = 0.0, tw = 0.0, wy = 0.0, wz = 0.0, w = 0.0;
        for (size_t q = 0; q < u.size(); ++q) {
            const size_t h = u[q];
            row[1 + od[h]] += pe[h];                                  // np.bincount(od[u], weights=pe[u]): member order
            tot += pe[h];
            tw += t[h] * pe[h];
            if (pe[h] >= m_par.hot_threshold * pe[u[best]]) {
                wy += m_pos[od[h]][1] * pe[h];
                wz += m_pos[od[h]][2] * pe[h];
                w += pe[h];
            }
            member.push_back((double) f);
            member.push_back((double) h);
        }
        double* s = &summary[f * 8];
        s[0] = (double) f;
        s[1] = tot;
        s[2] = w > 0 ? wy / w : 0.0;
        s[3] = w > 0 ? wz / w : 0.0;
        s[4] = (double) m_plane[od[u[0]]];
        s[5] = (double) u.size();
        s[6] = tot > 0 ? tw / tot : row[0];
        s[7] = (double) u[0];
    }
    auto tensors = std::make_shared<ITensor::vector>();
    tensors->push_back(make_tensor("opflash", flash, mcol));
    tensors->push_back(make_tensor("flash_summary", summary, 8));
    tensors->push_back(make_tensor("adj_membership", member, 2));
    tensors->push_back(ten);
    Configuration md = in->metadata();
    md["producer"] = "FdvdAdjOpFlashFinder";
    md["nflash"] = (Json::UInt64) nf;
    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), md, ITensor::shared_vector(tensors));
    log->debug("call={} ident={} hits={} flashes={}", m_count, in->ident(), nh, nf);
    ++m_count;
    return true;
}
