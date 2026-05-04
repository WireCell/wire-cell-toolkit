#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include <cmath>

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;

// init_tagger_info: reset a TaggerInfo struct to its default values.
// In the toolkit we rely on C++ default-member-initializers on the struct
// itself (see NeutrinoTaggerInfo.h), so value-initializing the struct is
// sufficient — no 1200-line assignment list needed.
void PatternAlgorithms::init_tagger_info(TaggerInfo& ti)
{
    ti = TaggerInfo{};
}

// fill_kine_tree: reconstruct per-particle and total neutrino kinetic energy
// starting from main_vertex, traversing the PR graph, and collecting all
// connected showers and track segments.
//
// Translation notes from prototype (WCPPID::NeutrinoID::fill_kine_tree):
//   map_vertex_segments[vtx]              -> boost::out_edges(vtx->get_descriptor(), graph) + graph[*ei].segment
//   find_other_vertex(seg, vtx)           -> find_other_vertex(graph, seg, vtx)
//   shower->get_start_segment()           -> shower->start_segment()
//   shower->get_start_segment()->get_particle_type() -> shower->get_particle_type()
//   shower->get_start_segment()->get_particle_mass() -> shower->start_segment()->particle_info()->mass()
//   shower->get_start_segment()->get_length()        -> segment_track_length(shower->start_segment())
//   seg->get_particle_type()              -> seg->particle_info()->pdg()
//   seg->get_kine_best()                  -> seg->particle_info()->kinetic_energy()
//   seg->get_particle_mass()              -> seg->particle_info()->mass()
//   cal_kine_charge(seg)                  -> cal_kine_charge(seg, graph, track_fitter, dv)
//   seg->cal_kine_dQdx()                  -> segment_cal_kine_dQdx(seg, recomb_model)
//   seg->cal_kine_range()                 -> cal_kine_range(segment_track_length(seg), pdg, particle_data)
//   shower->get_start_vertex()            -> shower->get_start_vertex_and_type()
//   kine_pio_*                            -> pio_kine struct fields
//   SCE correction                        -> geom_helper->get_corrected_point(...) (or raw point if null)
KineInfo PatternAlgorithms::fill_kine_tree(
    VertexPtr main_vertex,
    IndexedShowerSet& showers,
    const Pi0KineFeatures& pio_kine,
    Graph& graph,
    TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv,
    WireCell::IClusGeomHelper::pointer geom_helper,
    const Clus::ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model)
{
    KineInfo ktree{};

    // -------------------------------------------------------------------------
    // Neutrino vertex position with optional SCE correction
    // -------------------------------------------------------------------------
    Point nu_vtx = main_vertex->fit().point;

    if (geom_helper) {
        WirePlaneId wpid = dv->contained_by(nu_vtx);
        int apa  = wpid.apa();
        int face = wpid.face();
        Point corr = geom_helper->get_corrected_point(nu_vtx, IClusGeomHelper::SCE, apa, face);
        ktree.kine_nu_x_corr = static_cast<float>(corr.x() / units::cm);
        ktree.kine_nu_y_corr = static_cast<float>(corr.y() / units::cm);
        ktree.kine_nu_z_corr = static_cast<float>(corr.z() / units::cm);
    }
    else {
        // TODO: SCE correction requires a valid geom_helper; using raw vertex position for now.
        ktree.kine_nu_x_corr = static_cast<float>(nu_vtx.x() / units::cm);
        ktree.kine_nu_y_corr = static_cast<float>(nu_vtx.y() / units::cm);
        ktree.kine_nu_z_corr = static_cast<float>(nu_vtx.z() / units::cm);
    }

    // -------------------------------------------------------------------------
    // Build a map from a shower's start-segment to the shower itself, and
    // collect all vertices/segments already owned by showers.
    // -------------------------------------------------------------------------
    const double ave_binding_energy = 8.6 * units::MeV;

    IndexedVertexSet  used_vertices;
    IndexedSegmentSet used_segments;
    IndexedShowerSet  used_showers;

    // Mark all shower-internal vertices and segments as used.
    for (const ShowerPtr& shower : showers) {
        shower->fill_sets(used_vertices, used_segments, /*flag_exclude_start_segment=*/false);
    }

    // Map from start-segment -> shower, for fast lookup during graph traversal.
    std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp> map_sg_shower;
    for (const ShowerPtr& shower : showers) {
        SegmentPtr start_sg = shower->start_segment();
        map_sg_shower[start_sg] = shower;
    }

    // -------------------------------------------------------------------------
    // Helper: push one shower's kinematics into ktree vectors.
    // kine_energy_included is pushed by the caller (value differs by context).
    // -------------------------------------------------------------------------
    auto push_shower_kine = [&](const ShowerPtr& shower) {
        double kine_best   = shower->get_kine_best();
        double kine_charge = shower->get_kine_charge();
        double kine_range  = shower->get_kine_range();

        ktree.kine_energy_particle.push_back(static_cast<float>(kine_best / units::MeV));
        ktree.kine_particle_type.push_back(shower->get_particle_type());

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2); // charge
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1); // range
        else
            ktree.kine_energy_info.push_back(0); // dQdx

        // Add rest-mass correction for non-electrons/positrons
        if (shower->get_particle_type() != 11) {
            SegmentPtr start_sg = shower->start_segment();
            if (start_sg && start_sg->particle_info()) {
                ktree.kine_reco_add_energy += static_cast<float>(
                    start_sg->particle_info()->mass() / units::MeV);
            }
        }
    };

    // -------------------------------------------------------------------------
    // Helper: push one track segment's kinematics into ktree vectors.
    // Returns the segment's PDG code.
    // -------------------------------------------------------------------------
    auto push_segment_kine = [&](SegmentPtr seg, int include_flag) -> int {
        int    pdg        = 0;
        double mass       = 0;
        double kine_best  = 0;

        if (seg->particle_info()) {
            pdg       = seg->particle_info()->pdg();
            mass      = seg->particle_info()->mass();
            kine_best = seg->particle_info()->kinetic_energy();
        }
        double kine_charge = cal_kine_charge(seg, graph, track_fitter, dv);
        double kine_range  = cal_kine_range(segment_track_length(seg), pdg, particle_data);

        ktree.kine_energy_particle.push_back(static_cast<float>(kine_best / units::MeV));
        ktree.kine_particle_type.push_back(pdg);

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2);
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1);
        else
            ktree.kine_energy_info.push_back(0);

        ktree.kine_energy_included.push_back(include_flag);

        if (pdg == 2212) { // proton: add binding energy
            ktree.kine_reco_add_energy += static_cast<float>(ave_binding_energy / units::MeV);
        }
        else if (pdg != 11) { // not electron: add rest mass
            ktree.kine_reco_add_energy += static_cast<float>(mass / units::MeV);
        }
        return pdg;
    };

    // -------------------------------------------------------------------------
    // First pass: segments directly connected to main_vertex
    // -------------------------------------------------------------------------
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;

    auto [ei_begin, ei_end] = boost::out_edges(main_vertex->get_descriptor(), graph);
    for (auto ei = ei_begin; ei != ei_end; ++ei) {
        SegmentPtr seg = graph[*ei].segment;

        auto it = map_sg_shower.find(seg);
        if (it != map_sg_shower.end()) {
            // This segment is a shower start-segment.
            push_shower_kine(it->second);
            ktree.kine_energy_included.push_back(1);
            used_showers.insert(it->second);
        }
        else {
            // Track segment.
            used_segments.insert(seg);
            VertexPtr other_vtx = find_other_vertex(graph, seg, main_vertex);
            segments_to_be_examined.emplace_back(other_vtx, seg);
            push_segment_kine(seg, 1);
        }
    }
    used_vertices.insert(main_vertex);

    // -------------------------------------------------------------------------
    // BFS traversal of the remaining track graph
    // -------------------------------------------------------------------------
    while (!segments_to_be_examined.empty()) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;

        for (auto& [curr_vtx, prev_sg] : segments_to_be_examined) {
            if (used_vertices.count(curr_vtx)) continue;

            bool flag_reduce = false;
            int  prev_pdg = 0;
            if (prev_sg->particle_info()) prev_pdg = prev_sg->particle_info()->pdg();

            auto [ei2_begin, ei2_end] = boost::out_edges(curr_vtx->get_descriptor(), graph);
            for (auto ei2 = ei2_begin; ei2 != ei2_end; ++ei2) {
                SegmentPtr curr_sg = graph[*ei2].segment;
                if (curr_sg == prev_sg) continue;

                int curr_pdg = 0;
                if (curr_sg->particle_info()) curr_pdg = curr_sg->particle_info()->pdg();

                // Detect particle continuation (same type, or muon<->pion flip)
                if (curr_pdg == prev_pdg ||
                    (prev_pdg == 211 && curr_pdg == 13) ||
                    (prev_pdg == 13  && curr_pdg == 211))
                    flag_reduce = true;

                auto it2 = map_sg_shower.find(curr_sg);
                if (it2 == map_sg_shower.end()) {
                    // Track segment
                    if (used_segments.count(curr_sg)) continue;
                    used_segments.insert(curr_sg);

                    push_segment_kine(curr_sg, 1);

                    VertexPtr other_vtx = find_other_vertex(graph, curr_sg, curr_vtx);
                    if (!used_vertices.count(other_vtx))
                        temp_segments.emplace_back(other_vtx, curr_sg);
                }
                else {
                    // Shower
                    const ShowerPtr& shower = it2->second;
                    if (!used_showers.count(shower)) {
                        push_shower_kine(shower);
                        ktree.kine_energy_included.push_back(1);
                        used_showers.insert(shower);
                    }
                }
            }
            used_vertices.insert(curr_vtx);

            // If we detected a particle continuation, undo the rest-mass/binding-energy
            // added for prev_sg (it was already counted earlier in the chain).
            if (flag_reduce && prev_sg->particle_info()) {
                if (prev_pdg == 2212) {
                    ktree.kine_reco_add_energy -= static_cast<float>(ave_binding_energy / units::MeV);
                }
                else if (prev_pdg != 11) {
                    ktree.kine_reco_add_energy -= static_cast<float>(
                        prev_sg->particle_info()->mass() / units::MeV);
                }
            }
        }

        segments_to_be_examined = std::move(temp_segments);
    }

    // -------------------------------------------------------------------------
    // Remaining showers not yet attached to the traversal above
    // (e.g. secondary showers with start vertex type <= 3)
    // -------------------------------------------------------------------------
    for (const ShowerPtr& shower : showers) {
        if (used_showers.count(shower)) continue;

        auto [start_vtx, vtx_type] = shower->get_start_vertex_and_type();
        if (vtx_type > 3) continue;

        double kine_best   = shower->get_kine_best();
        double kine_charge = shower->get_kine_charge();
        double kine_range  = shower->get_kine_range();

        ktree.kine_energy_particle.push_back(static_cast<float>(kine_best / units::MeV));
        ktree.kine_particle_type.push_back(shower->get_particle_type());

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2);
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1);
        else
            ktree.kine_energy_info.push_back(0);

        ktree.kine_energy_included.push_back(vtx_type != 3 ? 1 : vtx_type);

        // Binding energy correction for proton showers with length > 5 cm
        if (shower->get_particle_type() == 2212) {
            SegmentPtr start_sg = shower->start_segment();
            if (start_sg && segment_track_length(start_sg) > 5.0 * units::cm) {
                ktree.kine_reco_add_energy += static_cast<float>(ave_binding_energy / units::MeV);
            }
        }
        // electrons: no rest-mass correction; other non-electrons not present in remaining showers
        // since push_shower_kine would have already handled them in the BFS phase

        used_showers.insert(shower);
    }

    // -------------------------------------------------------------------------
    // Total reconstructed neutrino energy
    // -------------------------------------------------------------------------
    ktree.kine_reco_Enu = 0;
    for (float e : ktree.kine_energy_particle) {
        ktree.kine_reco_Enu += e;
    }
    ktree.kine_reco_Enu += ktree.kine_reco_add_energy;

    // -------------------------------------------------------------------------
    // Pi0 kinematics (from pio_kine struct, angles converted to degrees)
    // -------------------------------------------------------------------------
    ktree.kine_pio_mass    = static_cast<float>(pio_kine.mass    / units::MeV);
    ktree.kine_pio_flag    = pio_kine.flag;
    ktree.kine_pio_vtx_dis = static_cast<float>(pio_kine.vtx_dis / units::cm);

    ktree.kine_pio_energy_1 = static_cast<float>(pio_kine.energy_1 / units::MeV);
    ktree.kine_pio_theta_1  = static_cast<float>(pio_kine.theta_1  / M_PI * 180.0);
    ktree.kine_pio_phi_1    = static_cast<float>(pio_kine.phi_1    / M_PI * 180.0);
    ktree.kine_pio_dis_1    = static_cast<float>(pio_kine.dis_1    / units::cm);

    ktree.kine_pio_energy_2 = static_cast<float>(pio_kine.energy_2 / units::MeV);
    ktree.kine_pio_theta_2  = static_cast<float>(pio_kine.theta_2  / M_PI * 180.0);
    ktree.kine_pio_phi_2    = static_cast<float>(pio_kine.phi_2    / M_PI * 180.0);
    ktree.kine_pio_dis_2    = static_cast<float>(pio_kine.dis_2    / units::cm);

    ktree.kine_pio_angle = static_cast<float>(pio_kine.angle / M_PI * 180.0);

    return ktree;
}
