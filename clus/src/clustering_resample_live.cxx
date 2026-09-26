// doc sbnd_xin/pr/149 round 2 -- ClusteringResampleLive.
//
// Re-sample every live blob's "3d" point cloud inside the pattern-recognition
// job, with the configured IBlobSampler, keeping clusters and blobs as they are.
//
// prototype (wire-cell-prod-nue.cxx:1289-1294, same step wire-cell-prod-stm.cxx:734):
//     // replace by the new sampling points ...
//     for (size_t i=0; i!=live_clusters.size();i++){
//       WCPPID::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
//       live_clusters.at(i)->Create_point_cloud();
//     }
// The prototype's clustering executable samples with the stepped rule
// (WCP2dToy::calc_sampling_points, 2dtoy/src/CalcPoints.cxx:613-700); its PR
// executables then replace every live cluster's cloud with the charge_stepped rule,
// disable_mix_dead_cell = true (pid/inc/WCPPID/CalcPoints.h:9-10), BEFORE
// Protect_Over_Clustering and NeutrinoID.  The toolkit's PR job instead reads the
// clustering job's stepped cloud from the saved tree.  This visitor is that missing
// step, as a pipeline stage; it must run FIRST (before switch_scope, which rebuilds
// x_t0cor and the in-volume split from the "3d" points).
//
// The PR job holds no IBlob / ISlice, so each is rebuilt from the tree:
//   shape    -- the blob's "scalar" wire bounds (ResampleLive::shape_from_bounds);
//   activity -- the grouping's ctpc row of the blob's slice (live), else the
//               dead-wind registry (dead), else absent (ResampleLive::compose_activity).
// With a 'stepped' sampler configured identically to the clustering job's, the
// rebuilt cloud must equal the saved one bit for bit: that is the identity gate of
// doc pr/149 round 2.
//
// Children are replaced as NODES, in their original order, not edited in place:
// Facade::Blob caches center/npoints (Facade_Blob.cxx fill_cache) and scoped views
// snapshot arrays; node removal/insertion is the path that notifies both.  The
// cluster-level "perblob" provenance rows stay aligned because the order is kept.
//
// Not in a pipeline => never constructed => the job is byte-identical.

#include "WireCellClus/ResampleLive.h"

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/WirePlaneId.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

class ClusteringResampleLive;
WIRECELL_FACTORY(ClusteringResampleLive, ClusteringResampleLive,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

namespace RL = WireCell::Clus::ResampleLive;

static Log::logptr_t logger()
{
    static Log::logptr_t l = Log::logger("clus.ResampleLive");
    return l;
}

class ClusteringResampleLive : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV {
  public:
    ClusteringResampleLive() {}
    virtual ~ClusteringResampleLive() {}

    virtual Configuration default_configuration() const
    {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["wire_margin"] = m_wire_margin;
        return cfg;
    }

    void configure(const WireCell::Configuration& cfg)
    {
        NeedDV::configure(cfg);
        m_grouping_name = get<std::string>(cfg, "grouping", m_grouping_name);
        m_wire_margin = get<int>(cfg, "wire_margin", m_wire_margin);

        for (const auto& aname : cfg["anodes"]) {
            auto anode = Factory::find_tn<IAnodePlane>(aname.asString());
            m_anode[anode->ident()] = anode;
            for (const auto& iface : anode->faces()) {
                if (!iface) continue;
                if (iface->raygrid().nlayers() != 5) {
                    raise<ValueError>("ClusteringResampleLive: unexpected number of ray grid layers %d",
                                      iface->raygrid().nlayers());
                }
                m_face[anode->ident()][iface->which()] = iface;
            }
        }
        // Same schema as RetileCluster (retile_cluster.cxx configure).
        for (const auto& scfg : cfg["samplers"]) {
            const int apa = scfg["apa"].asInt();
            const int face = scfg["face"].asInt();
            const std::string name = scfg["name"].asString();
            if (name.empty()) {
                raise<ValueError>("ClusteringResampleLive requires an IBlobSampler name for APA %d face %d", apa, face);
            }
            m_samplers[apa][face] = Factory::find_tn<IBlobSampler>(name);
        }
        if (m_samplers.empty()) {
            raise<ValueError>("ClusteringResampleLive: no samplers configured");
        }
    }

    void visit(Ensemble& ensemble) const;

  private:
    std::string m_grouping_name{"live"};
    // Activity is rebuilt for each plane's bounds widened by this many wires on
    // each side: the per-point wire-index extras and the third-plane round() can
    // land one wire past a strip bound.
    int m_wire_margin{2};
    std::map<int, IAnodePlane::pointer> m_anode;
    std::map<int, std::map<int, IAnodeFace::pointer>> m_face;
    std::map<int, std::map<int, IBlobSampler::pointer>> m_samplers;
};

namespace {
    struct Census {
        size_t clusters{0}, blobs{0}, points_before{0}, points_after{0};
        size_t kept_empty{0};      // resample gave no point: node copied unchanged
        size_t no_scalar{0};       // blob without a "scalar" PC: node copied unchanged
        std::array<size_t, 3> dead_edge{0, 0, 0};  // blobs with a dead first/last strip wire, per plane
        size_t dead_edge_gap_only{0};  // ... of which W only through the hand-declared gap registry
        size_t interval_mismatch{0};   // new max/min wire interval or type != saved
    };

    int scalar_int(const PointCloud::Dataset& ds, const char* name)
    {
        auto arr = ds.get(name);
        if (!arr) {
            raise<ValueError>("ClusteringResampleLive: blob scalar lacks '%s'", name);
        }
        return arr->elements<int>()[0];
    }
}

void ClusteringResampleLive::visit(Ensemble& ensemble) const
{
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        logger()->warn("ClusteringResampleLive: no '{}' grouping, skipping", m_grouping_name);
        return;
    }
    Grouping& grouping = *groupings.at(0);
    const auto ticks = grouping.get_tick();

    // Per-(apa, face) 2-D projection angles, computed as PointTreeBuilding does for
    // the clustering job's cloud (PointTreeBuilding.cxx operator(), angles loop).
    std::map<std::pair<int, int>, std::vector<double>> angles_by_face;
    auto angles_for = [&](int apa, int face) -> const std::vector<double>& {
        auto key = std::make_pair(apa, face);
        auto it = angles_by_face.find(key);
        if (it != angles_by_face.end()) return it->second;
        std::vector<double> angles(3);
        for (size_t ind = 0; ind < 3; ++ind) {
            WirePlaneId wpid(iplane2layer[ind], face, apa);
            Vector wire_dir = m_dv->wire_direction(wpid);
            angles[ind] = std::atan2(wire_dir.z(), wire_dir.y());
        }
        return angles_by_face.emplace(key, angles).first->second;
    };

    Census census;
    int blob_counter = 0;

    // Snapshot: nothing below adds or removes clusters.
    const std::vector<Cluster*> clusters = grouping.children();
    for (Cluster* cluster : clusters) {
        auto* cnode = cluster->node();
        ++census.clusters;

        // Detach every blob node first (emits the "removing" notices while the
        // nodes still hold their PCs, so scoped views containing them are erased),
        // then build replacements from the detached nodes, then insert in order.
        auto kids = cnode->remove_children();

        std::vector<named_pointclouds_t> replacements;
        replacements.reserve(kids.size());
        for (auto& kid : kids) {
            ++census.blobs;
            auto& lpcs = kid->value.local_pcs();
            auto sit = lpcs.find("scalar");
            auto pit = lpcs.find("3d");
            if (sit == lpcs.end()) {
                ++census.no_scalar;
                replacements.push_back(std::move(lpcs));
                continue;
            }
            const auto& scalar = sit->second;
            census.points_before += (pit == lpcs.end()) ? 0 : pit->second.size_major();

            const WirePlaneId wpid(scalar_int(scalar, "wpid"));
            const int apa = wpid.apa();
            const int face = wpid.face();
            const int smin = scalar_int(scalar, "slice_index_min");
            const int smax = scalar_int(scalar, "slice_index_max");
            const RL::plane_bounds_t bounds = {
                std::make_pair(scalar_int(scalar, "u_wire_index_min"), scalar_int(scalar, "u_wire_index_max")),
                std::make_pair(scalar_int(scalar, "v_wire_index_min"), scalar_int(scalar, "v_wire_index_max")),
                std::make_pair(scalar_int(scalar, "w_wire_index_min"), scalar_int(scalar, "w_wire_index_max")),
            };

            auto aface = m_face.find(apa);
            auto asamp = m_samplers.find(apa);
            if (aface == m_face.end() || !aface->second.count(face) ||
                asamp == m_samplers.end() || !asamp->second.count(face) || !m_anode.count(apa)) {
                raise<ValueError>("ClusteringResampleLive: no anode face / sampler for APA %d face %d", apa, face);
            }
            const auto& iface = aface->second.at(face);
            const auto& sampler = asamp->second.at(face);
            const auto& ianode = m_anode.at(apa);
            const double tick = ticks.at(apa).at(face);

            // Slice activity over each plane's bounds +- margin.
            ISlice::map_t activity;
            const auto& gap_w = grouping.get_dead_gap_winds(apa, face);
            bool any_dead_edge = false, gap_only = true;
            for (int p = 0; p < 3; ++p) {
                const auto& wires = iface->planes()[p]->wires();
                const int nwires = (int) wires.size();
                const int lo = std::max(0, bounds[p].first - m_wire_margin);
                const int hi = std::min(nwires - 1, bounds[p].second + m_wire_margin);
                const auto* row = grouping.wire_charge_row(apa, face, p, smin);
                auto live = [&](int w, double& q, double& err) {
                    if (!row) return false;
                    auto it = row->find(w);
                    if (it == row->end()) return false;
                    q = it->second.first;
                    err = it->second.second;
                    return true;
                };
                auto dead = [&](int w) { return grouping.is_wire_dead(apa, face, p, w, smin); };
                bool edge_dead = false, edge_gap = true;
                for (const auto& wv : RL::compose_activity(lo, hi, live, dead)) {
                    auto ich = ianode->channel(wires[wv.wire]->channel());
                    if (!ich) continue;
                    activity[ich] = ISlice::value_t(wv.charge, wv.error);
                    const bool is_edge = (wv.wire == bounds[p].first || wv.wire == bounds[p].second - 1);
                    if (is_edge && wv.error > 1e10) {   // BlobSampler dead_threshold default
                        edge_dead = true;
                        if (!(p == 2 && gap_w.count(wv.wire))) edge_gap = false;
                    }
                }
                if (edge_dead) {
                    ++census.dead_edge[p];
                    any_dead_edge = true;
                    if (!edge_gap) gap_only = false;
                }
            }
            if (any_dead_edge && gap_only) ++census.dead_edge_gap_only;

            auto islice = std::make_shared<Aux::SimpleSlice>(nullptr, smin, smin * tick, (smax - smin) * tick, activity);
            RayGrid::Blob shape = RL::shape_from_bounds(iface->raygrid(), bounds);
            auto iblob = std::make_shared<Aux::SimpleBlob>(blob_counter, 0.0f, 0.0f, shape, islice, iface);
            auto fresh = Aux::sample_live(sampler, iblob, angles_for(apa, face), tick, blob_counter);
            ++blob_counter;

            auto fit = fresh.find("3d");
            if (fit == fresh.end() || fit->second.size_major() == 0) {
                ++census.kept_empty;
                census.points_after += (pit == lpcs.end()) ? 0 : pit->second.size_major();
                replacements.push_back(std::move(lpcs));
                continue;
            }

            const auto& fscalar = fresh.at("scalar");
            for (const char* name : {"max_wire_interval", "min_wire_interval", "max_wire_type", "min_wire_type"}) {
                auto a = scalar.get(name);
                auto b = fscalar.get(name);
                if (a && b && a->elements<int>()[0] != b->elements<int>()[0]) {
                    ++census.interval_mismatch;
                    break;
                }
            }

            named_pointclouds_t out;
            for (auto& [name, pc] : lpcs) {
                if (name == "3d") continue;
                out.emplace(name, std::move(pc));
            }
            RL::refresh_scalar(out.at("scalar"), fscalar);
            census.points_after += fit->second.size_major();
            out.emplace("3d", std::move(fit->second));
            replacements.push_back(std::move(out));
        }
        kids.clear();   // destroy the old nodes (and their facades) before re-inserting

        for (auto& pcs : replacements) {
            cnode->insert(Points(std::move(pcs)));
        }
        cluster->invalidate_cache();
    }

    SPDLOG_LOGGER_DEBUG(logger(),
                        "RESAMPLE clusters {} blobs {} points {} -> {} kept_empty {} no_scalar {} "
                        "dead_edge_u {} dead_edge_v {} dead_edge_w {} dead_edge_gap_only {} interval_mismatch {}",
                        census.clusters, census.blobs, census.points_before, census.points_after,
                        census.kept_empty, census.no_scalar, census.dead_edge[0], census.dead_edge[1],
                        census.dead_edge[2], census.dead_edge_gap_only, census.interval_mismatch);
}
