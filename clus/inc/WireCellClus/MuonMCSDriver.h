/** MCS muon momentum driver (doc 80 round 2).
 *
 * Bridges the PR data model onto the ROOT-free WireCellMcs library: selects
 * the muon, harvests its trajectory cloud, converts internal units -> cm on
 * the boundary, runs Mcs::MuonMCS, and stamps the result into PR::KineInfo's
 * kine_mcs_* scalars.  Called from TaggerCheckNeutrino::visit() immediately
 * before track_fitter->set_kine_info(), once per bundle -- a knob-gated call
 * site, deliberately NOT a new visitor (doc 80 sec 7: a visitor could not see
 * the function-local segments_in_long_muon, and a pipeline_names entry would
 * move the compiled config even when off).
 */
#ifndef WIRECELLCLUS_MUONMCSDRIVER
#define WIRECELLCLUS_MUONMCSDRIVER

#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellUtil/Logging.h"

#include <string>

namespace WireCell::Clus::PR {

    struct MuonMCSConfig {
        // All defaults are the legacy no-op / owner-decided values of doc 80
        // sec 7.1.  seg_length (14 cm = X0(LAr)) is deliberately NOT here: it
        // is a structural constant of the tune, not a knob.
        bool enable{false};
        std::string muon_source{"pf_muon"};          // pf_muon | long_muon | longest_segment | long_muon_else_pf (doc 84 round 1)
        double muon_min_length_cm{40};
        std::string point_source{"muon_segments"};   // muon_segments | whole_event (validation only)
        bool beam_window_only{true};   // correctness, not cost: out-of-spill readout
                                       // may truncate the track, corrupting range AND
                                       // MCS silently (doc 80 sec 7.4)
        double cathode_x_cm{0};        // cathode plane x [cm]
        double cathode_xcut_cm{0};     // excised half-band [cm]; 0 = off (SBND: 5)
        int max_points{20000};         // whole_event perf guard
        bool range_comparator_chain{false};  // doc 84 round 1 (P5): extra log-only
                                             // sentinel with the long-muon chain's
                                             // summed-length range KE; no output
                                             // bytes move either way
    };

    /// Run MCS for this bundle and fill kine's kine_mcs_* fields.  Fields are
    /// left at their -1 defaults when no muon qualifies.  beam_gate_active:
    /// whether the candidate list was built under a beam-window gate (the
    /// per-bundle candidates are already spill-coincident in that case).
    /// particle_data/recomb_model (both may be null): when present, the INFO
    /// line additionally carries the toolkit's cal_kine_range and
    /// segment_cal_kine_dQdx KE for the SAME selected muon -- the round-4
    /// comparators, emitted as a log sentinel (pr/94 ROW precedent) because
    /// the T_kine schema deliberately stays at the five kine_mcs_* scalars.
    void mcs_fill_kine(KineInfo& kine, Graph& graph,
                       const IndexedSegmentSet& segments_in_long_muon,
                       const VertexPtr& main_vertex,
                       bool beam_gate_active,
                       const MuonMCSConfig& cfg,
                       const Clus::ParticleDataSet::pointer& particle_data,
                       const IRecombinationModel::pointer& recomb_model,
                       WireCell::Log::logptr_t log);

}  // namespace WireCell::Clus::PR

#endif
