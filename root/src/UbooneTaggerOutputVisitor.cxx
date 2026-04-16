#include "WireCellRoot/UbooneTaggerOutputVisitor.h"

#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"

WIRECELL_FACTORY(UbooneTaggerOutputVisitor, WireCell::Root::UbooneTaggerOutputVisitor,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;

Root::UbooneTaggerOutputVisitor::UbooneTaggerOutputVisitor()
  : log(Log::logger("tagger_output"))
{
}

Root::UbooneTaggerOutputVisitor::~UbooneTaggerOutputVisitor() {}

void Root::UbooneTaggerOutputVisitor::configure(const WireCell::Configuration& cfg)
{
    m_output_filename = get<std::string>(cfg, "output_filename", "tracking_proj.root");
    m_grouping_name = get<std::string>(cfg, "grouping", "live");
}

WireCell::Configuration Root::UbooneTaggerOutputVisitor::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = "tracking_proj.root";
    cfg["grouping"] = "live";
    return cfg;
}

void Root::UbooneTaggerOutputVisitor::visit(Clus::Facade::Ensemble& ensemble) const
{
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("UbooneTaggerOutputVisitor: no grouping '{}'", m_grouping_name);
        return;
    }

    auto& grouping = *groupings.at(0);
    auto tf = grouping.get_track_fitting();
    if (!tf) {
        log->debug("UbooneTaggerOutputVisitor: no TrackFitting in grouping");
        return;
    }

    // Make mutable copies for branch addresses (visit is const).
    PR::TaggerInfo ti = tf->get_tagger_info();
    PR::KineInfo ki = tf->get_kine_info();

    // Open existing ROOT file in UPDATE mode to add trees.
    TFile* output_tf = TFile::Open(m_output_filename.c_str(), "UPDATE");
    if (!output_tf || output_tf->IsZombie()) {
        log->error("UbooneTaggerOutputVisitor: cannot open {} for update", m_output_filename);
        return;
    }

    // ================================================================
    //  T_tagger tree — all TaggerInfo fields + vertex position
    // ================================================================
    TTree* t_tagger = new TTree("T_tagger", "T_tagger");
    t_tagger->SetDirectory(output_tf);

    // Vertex position from KineInfo (matches prototype's nu_x/nu_y/nu_z branches)
    t_tagger->Branch("nu_x", &ki.kine_nu_x_corr, "nu_x/F");
    t_tagger->Branch("nu_y", &ki.kine_nu_y_corr, "nu_y/F");
    t_tagger->Branch("nu_z", &ki.kine_nu_z_corr, "nu_z/F");

    // ---- cosmic tagger (top-level flag) ----
    t_tagger->Branch("cosmic_flag", &ti.cosmic_flag, "cosmic_flag/F");
    t_tagger->Branch("cosmic_n_solid_tracks", &ti.cosmic_n_solid_tracks, "cosmic_n_solid_tracks/F");
    t_tagger->Branch("cosmic_energy_main_showers", &ti.cosmic_energy_main_showers, "cosmic_energy_main_showers/F");
    t_tagger->Branch("cosmic_energy_direct_showers", &ti.cosmic_energy_direct_showers, "cosmic_energy_direct_showers/F");
    t_tagger->Branch("cosmic_energy_indirect_showers", &ti.cosmic_energy_indirect_showers, "cosmic_energy_indirect_showers/F");
    t_tagger->Branch("cosmic_n_direct_showers", &ti.cosmic_n_direct_showers, "cosmic_n_direct_showers/F");
    t_tagger->Branch("cosmic_n_indirect_showers", &ti.cosmic_n_indirect_showers, "cosmic_n_indirect_showers/F");
    t_tagger->Branch("cosmic_n_main_showers", &ti.cosmic_n_main_showers, "cosmic_n_main_showers/F");
    t_tagger->Branch("cosmic_filled", &ti.cosmic_filled, "cosmic_filled/F");

    // ---- shower gap identification ----
    t_tagger->Branch("gap_flag", &ti.gap_flag, "gap_flag/F");
    t_tagger->Branch("gap_flag_prolong_u", &ti.gap_flag_prolong_u, "gap_flag_prolong_u/F");
    t_tagger->Branch("gap_flag_prolong_v", &ti.gap_flag_prolong_v, "gap_flag_prolong_v/F");
    t_tagger->Branch("gap_flag_prolong_w", &ti.gap_flag_prolong_w, "gap_flag_prolong_w/F");
    t_tagger->Branch("gap_flag_parallel", &ti.gap_flag_parallel, "gap_flag_parallel/F");
    t_tagger->Branch("gap_n_points", &ti.gap_n_points, "gap_n_points/F");
    t_tagger->Branch("gap_n_bad", &ti.gap_n_bad, "gap_n_bad/F");
    t_tagger->Branch("gap_energy", &ti.gap_energy, "gap_energy/F");
    t_tagger->Branch("gap_num_valid_tracks", &ti.gap_num_valid_tracks, "gap_num_valid_tracks/F");
    t_tagger->Branch("gap_flag_single_shower", &ti.gap_flag_single_shower, "gap_flag_single_shower/F");
    t_tagger->Branch("gap_filled", &ti.gap_filled, "gap_filled/F");

    // ---- MIP quality ----
    t_tagger->Branch("mip_quality_flag", &ti.mip_quality_flag, "mip_quality_flag/F");
    t_tagger->Branch("mip_quality_energy", &ti.mip_quality_energy, "mip_quality_energy/F");
    t_tagger->Branch("mip_quality_overlap", &ti.mip_quality_overlap, "mip_quality_overlap/F");
    t_tagger->Branch("mip_quality_n_showers", &ti.mip_quality_n_showers, "mip_quality_n_showers/F");
    t_tagger->Branch("mip_quality_n_tracks", &ti.mip_quality_n_tracks, "mip_quality_n_tracks/F");
    t_tagger->Branch("mip_quality_flag_inside_pi0", &ti.mip_quality_flag_inside_pi0, "mip_quality_flag_inside_pi0/F");
    t_tagger->Branch("mip_quality_n_pi0_showers", &ti.mip_quality_n_pi0_showers, "mip_quality_n_pi0_showers/F");
    t_tagger->Branch("mip_quality_shortest_length", &ti.mip_quality_shortest_length, "mip_quality_shortest_length/F");
    t_tagger->Branch("mip_quality_acc_length", &ti.mip_quality_acc_length, "mip_quality_acc_length/F");
    t_tagger->Branch("mip_quality_shortest_angle", &ti.mip_quality_shortest_angle, "mip_quality_shortest_angle/F");
    t_tagger->Branch("mip_quality_flag_proton", &ti.mip_quality_flag_proton, "mip_quality_flag_proton/F");
    t_tagger->Branch("mip_quality_filled", &ti.mip_quality_filled, "mip_quality_filled/F");

    // ---- MIP identification ----
    t_tagger->Branch("mip_flag", &ti.mip_flag, "mip_flag/F");
    t_tagger->Branch("mip_energy", &ti.mip_energy, "mip_energy/F");
    t_tagger->Branch("mip_n_end_reduction", &ti.mip_n_end_reduction, "mip_n_end_reduction/F");
    t_tagger->Branch("mip_n_first_mip", &ti.mip_n_first_mip, "mip_n_first_mip/F");
    t_tagger->Branch("mip_n_first_non_mip", &ti.mip_n_first_non_mip, "mip_n_first_non_mip/F");
    t_tagger->Branch("mip_n_first_non_mip_1", &ti.mip_n_first_non_mip_1, "mip_n_first_non_mip_1/F");
    t_tagger->Branch("mip_n_first_non_mip_2", &ti.mip_n_first_non_mip_2, "mip_n_first_non_mip_2/F");
    t_tagger->Branch("mip_vec_dQ_dx_0", &ti.mip_vec_dQ_dx_0, "mip_vec_dQ_dx_0/F");
    t_tagger->Branch("mip_vec_dQ_dx_1", &ti.mip_vec_dQ_dx_1, "mip_vec_dQ_dx_1/F");
    t_tagger->Branch("mip_vec_dQ_dx_2", &ti.mip_vec_dQ_dx_2, "mip_vec_dQ_dx_2/F");
    t_tagger->Branch("mip_vec_dQ_dx_3", &ti.mip_vec_dQ_dx_3, "mip_vec_dQ_dx_3/F");
    t_tagger->Branch("mip_vec_dQ_dx_4", &ti.mip_vec_dQ_dx_4, "mip_vec_dQ_dx_4/F");
    t_tagger->Branch("mip_vec_dQ_dx_5", &ti.mip_vec_dQ_dx_5, "mip_vec_dQ_dx_5/F");
    t_tagger->Branch("mip_vec_dQ_dx_6", &ti.mip_vec_dQ_dx_6, "mip_vec_dQ_dx_6/F");
    t_tagger->Branch("mip_vec_dQ_dx_7", &ti.mip_vec_dQ_dx_7, "mip_vec_dQ_dx_7/F");
    t_tagger->Branch("mip_vec_dQ_dx_8", &ti.mip_vec_dQ_dx_8, "mip_vec_dQ_dx_8/F");
    t_tagger->Branch("mip_vec_dQ_dx_9", &ti.mip_vec_dQ_dx_9, "mip_vec_dQ_dx_9/F");
    t_tagger->Branch("mip_vec_dQ_dx_10", &ti.mip_vec_dQ_dx_10, "mip_vec_dQ_dx_10/F");
    t_tagger->Branch("mip_vec_dQ_dx_11", &ti.mip_vec_dQ_dx_11, "mip_vec_dQ_dx_11/F");
    t_tagger->Branch("mip_vec_dQ_dx_12", &ti.mip_vec_dQ_dx_12, "mip_vec_dQ_dx_12/F");
    t_tagger->Branch("mip_vec_dQ_dx_13", &ti.mip_vec_dQ_dx_13, "mip_vec_dQ_dx_13/F");
    t_tagger->Branch("mip_vec_dQ_dx_14", &ti.mip_vec_dQ_dx_14, "mip_vec_dQ_dx_14/F");
    t_tagger->Branch("mip_vec_dQ_dx_15", &ti.mip_vec_dQ_dx_15, "mip_vec_dQ_dx_15/F");
    t_tagger->Branch("mip_vec_dQ_dx_16", &ti.mip_vec_dQ_dx_16, "mip_vec_dQ_dx_16/F");
    t_tagger->Branch("mip_vec_dQ_dx_17", &ti.mip_vec_dQ_dx_17, "mip_vec_dQ_dx_17/F");
    t_tagger->Branch("mip_vec_dQ_dx_18", &ti.mip_vec_dQ_dx_18, "mip_vec_dQ_dx_18/F");
    t_tagger->Branch("mip_vec_dQ_dx_19", &ti.mip_vec_dQ_dx_19, "mip_vec_dQ_dx_19/F");
    t_tagger->Branch("mip_max_dQ_dx_sample", &ti.mip_max_dQ_dx_sample, "mip_max_dQ_dx_sample/F");
    t_tagger->Branch("mip_n_below_threshold", &ti.mip_n_below_threshold, "mip_n_below_threshold/F");
    t_tagger->Branch("mip_n_below_zero", &ti.mip_n_below_zero, "mip_n_below_zero/F");
    t_tagger->Branch("mip_n_lowest", &ti.mip_n_lowest, "mip_n_lowest/F");
    t_tagger->Branch("mip_n_highest", &ti.mip_n_highest, "mip_n_highest/F");
    t_tagger->Branch("mip_lowest_dQ_dx", &ti.mip_lowest_dQ_dx, "mip_lowest_dQ_dx/F");
    t_tagger->Branch("mip_highest_dQ_dx", &ti.mip_highest_dQ_dx, "mip_highest_dQ_dx/F");
    t_tagger->Branch("mip_medium_dQ_dx", &ti.mip_medium_dQ_dx, "mip_medium_dQ_dx/F");
    t_tagger->Branch("mip_stem_length", &ti.mip_stem_length, "mip_stem_length/F");
    t_tagger->Branch("mip_length_main", &ti.mip_length_main, "mip_length_main/F");
    t_tagger->Branch("mip_length_total", &ti.mip_length_total, "mip_length_total/F");
    t_tagger->Branch("mip_angle_beam", &ti.mip_angle_beam, "mip_angle_beam/F");
    t_tagger->Branch("mip_iso_angle", &ti.mip_iso_angle, "mip_iso_angle/F");
    t_tagger->Branch("mip_n_vertex", &ti.mip_n_vertex, "mip_n_vertex/F");
    t_tagger->Branch("mip_n_good_tracks", &ti.mip_n_good_tracks, "mip_n_good_tracks/F");
    t_tagger->Branch("mip_E_indirect_max_energy", &ti.mip_E_indirect_max_energy, "mip_E_indirect_max_energy/F");
    t_tagger->Branch("mip_flag_all_above", &ti.mip_flag_all_above, "mip_flag_all_above/F");
    t_tagger->Branch("mip_min_dQ_dx_5", &ti.mip_min_dQ_dx_5, "mip_min_dQ_dx_5/F");
    t_tagger->Branch("mip_n_other_vertex", &ti.mip_n_other_vertex, "mip_n_other_vertex/F");
    t_tagger->Branch("mip_n_stem_size", &ti.mip_n_stem_size, "mip_n_stem_size/F");
    t_tagger->Branch("mip_flag_stem_trajectory", &ti.mip_flag_stem_trajectory, "mip_flag_stem_trajectory/F");
    t_tagger->Branch("mip_min_dis", &ti.mip_min_dis, "mip_min_dis/F");
    t_tagger->Branch("mip_filled", &ti.mip_filled, "mip_filled/F");

    // ---- SSM (short straight muon / KDAR) tagger ----
    t_tagger->Branch("ssm_flag_st_kdar", &ti.ssm_flag_st_kdar, "ssm_flag_st_kdar/F");
    t_tagger->Branch("ssm_Nsm", &ti.ssm_Nsm, "ssm_Nsm/F");
    t_tagger->Branch("ssm_Nsm_wivtx", &ti.ssm_Nsm_wivtx, "ssm_Nsm_wivtx/F");
    t_tagger->Branch("ssm_dq_dx_fwd_1", &ti.ssm_dq_dx_fwd_1, "ssm_dq_dx_fwd_1/F");
    t_tagger->Branch("ssm_dq_dx_fwd_2", &ti.ssm_dq_dx_fwd_2, "ssm_dq_dx_fwd_2/F");
    t_tagger->Branch("ssm_dq_dx_fwd_3", &ti.ssm_dq_dx_fwd_3, "ssm_dq_dx_fwd_3/F");
    t_tagger->Branch("ssm_dq_dx_fwd_4", &ti.ssm_dq_dx_fwd_4, "ssm_dq_dx_fwd_4/F");
    t_tagger->Branch("ssm_dq_dx_fwd_5", &ti.ssm_dq_dx_fwd_5, "ssm_dq_dx_fwd_5/F");
    t_tagger->Branch("ssm_dq_dx_bck_1", &ti.ssm_dq_dx_bck_1, "ssm_dq_dx_bck_1/F");
    t_tagger->Branch("ssm_dq_dx_bck_2", &ti.ssm_dq_dx_bck_2, "ssm_dq_dx_bck_2/F");
    t_tagger->Branch("ssm_dq_dx_bck_3", &ti.ssm_dq_dx_bck_3, "ssm_dq_dx_bck_3/F");
    t_tagger->Branch("ssm_dq_dx_bck_4", &ti.ssm_dq_dx_bck_4, "ssm_dq_dx_bck_4/F");
    t_tagger->Branch("ssm_dq_dx_bck_5", &ti.ssm_dq_dx_bck_5, "ssm_dq_dx_bck_5/F");
    t_tagger->Branch("ssm_d_dq_dx_fwd_12", &ti.ssm_d_dq_dx_fwd_12, "ssm_d_dq_dx_fwd_12/F");
    t_tagger->Branch("ssm_d_dq_dx_fwd_23", &ti.ssm_d_dq_dx_fwd_23, "ssm_d_dq_dx_fwd_23/F");
    t_tagger->Branch("ssm_d_dq_dx_fwd_34", &ti.ssm_d_dq_dx_fwd_34, "ssm_d_dq_dx_fwd_34/F");
    t_tagger->Branch("ssm_d_dq_dx_fwd_45", &ti.ssm_d_dq_dx_fwd_45, "ssm_d_dq_dx_fwd_45/F");
    t_tagger->Branch("ssm_d_dq_dx_bck_12", &ti.ssm_d_dq_dx_bck_12, "ssm_d_dq_dx_bck_12/F");
    t_tagger->Branch("ssm_d_dq_dx_bck_23", &ti.ssm_d_dq_dx_bck_23, "ssm_d_dq_dx_bck_23/F");
    t_tagger->Branch("ssm_d_dq_dx_bck_34", &ti.ssm_d_dq_dx_bck_34, "ssm_d_dq_dx_bck_34/F");
    t_tagger->Branch("ssm_d_dq_dx_bck_45", &ti.ssm_d_dq_dx_bck_45, "ssm_d_dq_dx_bck_45/F");
    t_tagger->Branch("ssm_max_dq_dx_fwd_3", &ti.ssm_max_dq_dx_fwd_3, "ssm_max_dq_dx_fwd_3/F");
    t_tagger->Branch("ssm_max_dq_dx_fwd_5", &ti.ssm_max_dq_dx_fwd_5, "ssm_max_dq_dx_fwd_5/F");
    t_tagger->Branch("ssm_max_dq_dx_bck_3", &ti.ssm_max_dq_dx_bck_3, "ssm_max_dq_dx_bck_3/F");
    t_tagger->Branch("ssm_max_dq_dx_bck_5", &ti.ssm_max_dq_dx_bck_5, "ssm_max_dq_dx_bck_5/F");
    t_tagger->Branch("ssm_max_d_dq_dx_fwd_3", &ti.ssm_max_d_dq_dx_fwd_3, "ssm_max_d_dq_dx_fwd_3/F");
    t_tagger->Branch("ssm_max_d_dq_dx_fwd_5", &ti.ssm_max_d_dq_dx_fwd_5, "ssm_max_d_dq_dx_fwd_5/F");
    t_tagger->Branch("ssm_max_d_dq_dx_bck_3", &ti.ssm_max_d_dq_dx_bck_3, "ssm_max_d_dq_dx_bck_3/F");
    t_tagger->Branch("ssm_max_d_dq_dx_bck_5", &ti.ssm_max_d_dq_dx_bck_5, "ssm_max_d_dq_dx_bck_5/F");
    t_tagger->Branch("ssm_medium_dq_dx", &ti.ssm_medium_dq_dx, "ssm_medium_dq_dx/F");
    t_tagger->Branch("ssm_medium_dq_dx_bp", &ti.ssm_medium_dq_dx_bp, "ssm_medium_dq_dx_bp/F");
    t_tagger->Branch("ssm_angle_to_z", &ti.ssm_angle_to_z, "ssm_angle_to_z/F");
    t_tagger->Branch("ssm_angle_to_target", &ti.ssm_angle_to_target, "ssm_angle_to_target/F");
    t_tagger->Branch("ssm_angle_to_absorber", &ti.ssm_angle_to_absorber, "ssm_angle_to_absorber/F");
    t_tagger->Branch("ssm_angle_to_vertical", &ti.ssm_angle_to_vertical, "ssm_angle_to_vertical/F");
    t_tagger->Branch("ssm_x_dir", &ti.ssm_x_dir, "ssm_x_dir/F");
    t_tagger->Branch("ssm_y_dir", &ti.ssm_y_dir, "ssm_y_dir/F");
    t_tagger->Branch("ssm_z_dir", &ti.ssm_z_dir, "ssm_z_dir/F");
    t_tagger->Branch("ssm_kine_energy", &ti.ssm_kine_energy, "ssm_kine_energy/F");
    t_tagger->Branch("ssm_kine_energy_reduced", &ti.ssm_kine_energy_reduced, "ssm_kine_energy_reduced/F");
    t_tagger->Branch("ssm_vtx_activity", &ti.ssm_vtx_activity, "ssm_vtx_activity/F");
    t_tagger->Branch("ssm_pdg", &ti.ssm_pdg, "ssm_pdg/F");
    t_tagger->Branch("ssm_dQ_dx_cut", &ti.ssm_dQ_dx_cut, "ssm_dQ_dx_cut/F");
    t_tagger->Branch("ssm_score_mu_fwd", &ti.ssm_score_mu_fwd, "ssm_score_mu_fwd/F");
    t_tagger->Branch("ssm_score_p_fwd", &ti.ssm_score_p_fwd, "ssm_score_p_fwd/F");
    t_tagger->Branch("ssm_score_e_fwd", &ti.ssm_score_e_fwd, "ssm_score_e_fwd/F");
    t_tagger->Branch("ssm_score_mu_bck", &ti.ssm_score_mu_bck, "ssm_score_mu_bck/F");
    t_tagger->Branch("ssm_score_p_bck", &ti.ssm_score_p_bck, "ssm_score_p_bck/F");
    t_tagger->Branch("ssm_score_e_bck", &ti.ssm_score_e_bck, "ssm_score_e_bck/F");
    t_tagger->Branch("ssm_score_mu_fwd_bp", &ti.ssm_score_mu_fwd_bp, "ssm_score_mu_fwd_bp/F");
    t_tagger->Branch("ssm_score_p_fwd_bp", &ti.ssm_score_p_fwd_bp, "ssm_score_p_fwd_bp/F");
    t_tagger->Branch("ssm_score_e_fwd_bp", &ti.ssm_score_e_fwd_bp, "ssm_score_e_fwd_bp/F");
    t_tagger->Branch("ssm_length", &ti.ssm_length, "ssm_length/F");
    t_tagger->Branch("ssm_direct_length", &ti.ssm_direct_length, "ssm_direct_length/F");
    t_tagger->Branch("ssm_length_ratio", &ti.ssm_length_ratio, "ssm_length_ratio/F");
    t_tagger->Branch("ssm_max_dev", &ti.ssm_max_dev, "ssm_max_dev/F");
    t_tagger->Branch("ssm_n_prim_tracks_1", &ti.ssm_n_prim_tracks_1, "ssm_n_prim_tracks_1/F");
    t_tagger->Branch("ssm_n_prim_tracks_3", &ti.ssm_n_prim_tracks_3, "ssm_n_prim_tracks_3/F");
    t_tagger->Branch("ssm_n_prim_tracks_5", &ti.ssm_n_prim_tracks_5, "ssm_n_prim_tracks_5/F");
    t_tagger->Branch("ssm_n_prim_tracks_8", &ti.ssm_n_prim_tracks_8, "ssm_n_prim_tracks_8/F");
    t_tagger->Branch("ssm_n_prim_tracks_11", &ti.ssm_n_prim_tracks_11, "ssm_n_prim_tracks_11/F");
    t_tagger->Branch("ssm_n_all_tracks_1", &ti.ssm_n_all_tracks_1, "ssm_n_all_tracks_1/F");
    t_tagger->Branch("ssm_n_all_tracks_3", &ti.ssm_n_all_tracks_3, "ssm_n_all_tracks_3/F");
    t_tagger->Branch("ssm_n_all_tracks_5", &ti.ssm_n_all_tracks_5, "ssm_n_all_tracks_5/F");
    t_tagger->Branch("ssm_n_all_tracks_8", &ti.ssm_n_all_tracks_8, "ssm_n_all_tracks_8/F");
    t_tagger->Branch("ssm_n_all_tracks_11", &ti.ssm_n_all_tracks_11, "ssm_n_all_tracks_11/F");
    t_tagger->Branch("ssm_n_daughter_tracks_1", &ti.ssm_n_daughter_tracks_1, "ssm_n_daughter_tracks_1/F");
    t_tagger->Branch("ssm_n_daughter_tracks_3", &ti.ssm_n_daughter_tracks_3, "ssm_n_daughter_tracks_3/F");
    t_tagger->Branch("ssm_n_daughter_tracks_5", &ti.ssm_n_daughter_tracks_5, "ssm_n_daughter_tracks_5/F");
    t_tagger->Branch("ssm_n_daughter_tracks_8", &ti.ssm_n_daughter_tracks_8, "ssm_n_daughter_tracks_8/F");
    t_tagger->Branch("ssm_n_daughter_tracks_11", &ti.ssm_n_daughter_tracks_11, "ssm_n_daughter_tracks_11/F");
    t_tagger->Branch("ssm_n_daughter_all_1", &ti.ssm_n_daughter_all_1, "ssm_n_daughter_all_1/F");
    t_tagger->Branch("ssm_n_daughter_all_3", &ti.ssm_n_daughter_all_3, "ssm_n_daughter_all_3/F");
    t_tagger->Branch("ssm_n_daughter_all_5", &ti.ssm_n_daughter_all_5, "ssm_n_daughter_all_5/F");
    t_tagger->Branch("ssm_n_daughter_all_8", &ti.ssm_n_daughter_all_8, "ssm_n_daughter_all_8/F");
    t_tagger->Branch("ssm_n_daughter_all_11", &ti.ssm_n_daughter_all_11, "ssm_n_daughter_all_11/F");

    // SSM primary track 1
    t_tagger->Branch("ssm_prim_track1_pdg", &ti.ssm_prim_track1_pdg, "ssm_prim_track1_pdg/F");
    t_tagger->Branch("ssm_prim_track1_score_mu_fwd", &ti.ssm_prim_track1_score_mu_fwd, "ssm_prim_track1_score_mu_fwd/F");
    t_tagger->Branch("ssm_prim_track1_score_p_fwd", &ti.ssm_prim_track1_score_p_fwd, "ssm_prim_track1_score_p_fwd/F");
    t_tagger->Branch("ssm_prim_track1_score_e_fwd", &ti.ssm_prim_track1_score_e_fwd, "ssm_prim_track1_score_e_fwd/F");
    t_tagger->Branch("ssm_prim_track1_score_mu_bck", &ti.ssm_prim_track1_score_mu_bck, "ssm_prim_track1_score_mu_bck/F");
    t_tagger->Branch("ssm_prim_track1_score_p_bck", &ti.ssm_prim_track1_score_p_bck, "ssm_prim_track1_score_p_bck/F");
    t_tagger->Branch("ssm_prim_track1_score_e_bck", &ti.ssm_prim_track1_score_e_bck, "ssm_prim_track1_score_e_bck/F");
    t_tagger->Branch("ssm_prim_track1_length", &ti.ssm_prim_track1_length, "ssm_prim_track1_length/F");
    t_tagger->Branch("ssm_prim_track1_direct_length", &ti.ssm_prim_track1_direct_length, "ssm_prim_track1_direct_length/F");
    t_tagger->Branch("ssm_prim_track1_length_ratio", &ti.ssm_prim_track1_length_ratio, "ssm_prim_track1_length_ratio/F");
    t_tagger->Branch("ssm_prim_track1_max_dev", &ti.ssm_prim_track1_max_dev, "ssm_prim_track1_max_dev/F");
    t_tagger->Branch("ssm_prim_track1_kine_energy_range", &ti.ssm_prim_track1_kine_energy_range, "ssm_prim_track1_kine_energy_range/F");
    t_tagger->Branch("ssm_prim_track1_kine_energy_range_mu", &ti.ssm_prim_track1_kine_energy_range_mu, "ssm_prim_track1_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_prim_track1_kine_energy_range_p", &ti.ssm_prim_track1_kine_energy_range_p, "ssm_prim_track1_kine_energy_range_p/F");
    t_tagger->Branch("ssm_prim_track1_kine_energy_range_e", &ti.ssm_prim_track1_kine_energy_range_e, "ssm_prim_track1_kine_energy_range_e/F");
    t_tagger->Branch("ssm_prim_track1_kine_energy_cal", &ti.ssm_prim_track1_kine_energy_cal, "ssm_prim_track1_kine_energy_cal/F");
    t_tagger->Branch("ssm_prim_track1_medium_dq_dx", &ti.ssm_prim_track1_medium_dq_dx, "ssm_prim_track1_medium_dq_dx/F");
    t_tagger->Branch("ssm_prim_track1_x_dir", &ti.ssm_prim_track1_x_dir, "ssm_prim_track1_x_dir/F");
    t_tagger->Branch("ssm_prim_track1_y_dir", &ti.ssm_prim_track1_y_dir, "ssm_prim_track1_y_dir/F");
    t_tagger->Branch("ssm_prim_track1_z_dir", &ti.ssm_prim_track1_z_dir, "ssm_prim_track1_z_dir/F");
    t_tagger->Branch("ssm_prim_track1_add_daught_track_counts_1", &ti.ssm_prim_track1_add_daught_track_counts_1, "ssm_prim_track1_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_prim_track1_add_daught_all_counts_1", &ti.ssm_prim_track1_add_daught_all_counts_1, "ssm_prim_track1_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_prim_track1_add_daught_track_counts_5", &ti.ssm_prim_track1_add_daught_track_counts_5, "ssm_prim_track1_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_prim_track1_add_daught_all_counts_5", &ti.ssm_prim_track1_add_daught_all_counts_5, "ssm_prim_track1_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_prim_track1_add_daught_track_counts_11", &ti.ssm_prim_track1_add_daught_track_counts_11, "ssm_prim_track1_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_prim_track1_add_daught_all_counts_11", &ti.ssm_prim_track1_add_daught_all_counts_11, "ssm_prim_track1_add_daught_all_counts_11/F");

    // SSM primary track 2
    t_tagger->Branch("ssm_prim_track2_pdg", &ti.ssm_prim_track2_pdg, "ssm_prim_track2_pdg/F");
    t_tagger->Branch("ssm_prim_track2_score_mu_fwd", &ti.ssm_prim_track2_score_mu_fwd, "ssm_prim_track2_score_mu_fwd/F");
    t_tagger->Branch("ssm_prim_track2_score_p_fwd", &ti.ssm_prim_track2_score_p_fwd, "ssm_prim_track2_score_p_fwd/F");
    t_tagger->Branch("ssm_prim_track2_score_e_fwd", &ti.ssm_prim_track2_score_e_fwd, "ssm_prim_track2_score_e_fwd/F");
    t_tagger->Branch("ssm_prim_track2_score_mu_bck", &ti.ssm_prim_track2_score_mu_bck, "ssm_prim_track2_score_mu_bck/F");
    t_tagger->Branch("ssm_prim_track2_score_p_bck", &ti.ssm_prim_track2_score_p_bck, "ssm_prim_track2_score_p_bck/F");
    t_tagger->Branch("ssm_prim_track2_score_e_bck", &ti.ssm_prim_track2_score_e_bck, "ssm_prim_track2_score_e_bck/F");
    t_tagger->Branch("ssm_prim_track2_length", &ti.ssm_prim_track2_length, "ssm_prim_track2_length/F");
    t_tagger->Branch("ssm_prim_track2_direct_length", &ti.ssm_prim_track2_direct_length, "ssm_prim_track2_direct_length/F");
    t_tagger->Branch("ssm_prim_track2_length_ratio", &ti.ssm_prim_track2_length_ratio, "ssm_prim_track2_length_ratio/F");
    t_tagger->Branch("ssm_prim_track2_max_dev", &ti.ssm_prim_track2_max_dev, "ssm_prim_track2_max_dev/F");
    t_tagger->Branch("ssm_prim_track2_kine_energy_range", &ti.ssm_prim_track2_kine_energy_range, "ssm_prim_track2_kine_energy_range/F");
    t_tagger->Branch("ssm_prim_track2_kine_energy_range_mu", &ti.ssm_prim_track2_kine_energy_range_mu, "ssm_prim_track2_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_prim_track2_kine_energy_range_p", &ti.ssm_prim_track2_kine_energy_range_p, "ssm_prim_track2_kine_energy_range_p/F");
    t_tagger->Branch("ssm_prim_track2_kine_energy_range_e", &ti.ssm_prim_track2_kine_energy_range_e, "ssm_prim_track2_kine_energy_range_e/F");
    t_tagger->Branch("ssm_prim_track2_kine_energy_cal", &ti.ssm_prim_track2_kine_energy_cal, "ssm_prim_track2_kine_energy_cal/F");
    t_tagger->Branch("ssm_prim_track2_medium_dq_dx", &ti.ssm_prim_track2_medium_dq_dx, "ssm_prim_track2_medium_dq_dx/F");
    t_tagger->Branch("ssm_prim_track2_x_dir", &ti.ssm_prim_track2_x_dir, "ssm_prim_track2_x_dir/F");
    t_tagger->Branch("ssm_prim_track2_y_dir", &ti.ssm_prim_track2_y_dir, "ssm_prim_track2_y_dir/F");
    t_tagger->Branch("ssm_prim_track2_z_dir", &ti.ssm_prim_track2_z_dir, "ssm_prim_track2_z_dir/F");
    t_tagger->Branch("ssm_prim_track2_add_daught_track_counts_1", &ti.ssm_prim_track2_add_daught_track_counts_1, "ssm_prim_track2_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_prim_track2_add_daught_all_counts_1", &ti.ssm_prim_track2_add_daught_all_counts_1, "ssm_prim_track2_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_prim_track2_add_daught_track_counts_5", &ti.ssm_prim_track2_add_daught_track_counts_5, "ssm_prim_track2_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_prim_track2_add_daught_all_counts_5", &ti.ssm_prim_track2_add_daught_all_counts_5, "ssm_prim_track2_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_prim_track2_add_daught_track_counts_11", &ti.ssm_prim_track2_add_daught_track_counts_11, "ssm_prim_track2_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_prim_track2_add_daught_all_counts_11", &ti.ssm_prim_track2_add_daught_all_counts_11, "ssm_prim_track2_add_daught_all_counts_11/F");

    // SSM daughter track 1
    t_tagger->Branch("ssm_daught_track1_pdg", &ti.ssm_daught_track1_pdg, "ssm_daught_track1_pdg/F");
    t_tagger->Branch("ssm_daught_track1_score_mu_fwd", &ti.ssm_daught_track1_score_mu_fwd, "ssm_daught_track1_score_mu_fwd/F");
    t_tagger->Branch("ssm_daught_track1_score_p_fwd", &ti.ssm_daught_track1_score_p_fwd, "ssm_daught_track1_score_p_fwd/F");
    t_tagger->Branch("ssm_daught_track1_score_e_fwd", &ti.ssm_daught_track1_score_e_fwd, "ssm_daught_track1_score_e_fwd/F");
    t_tagger->Branch("ssm_daught_track1_score_mu_bck", &ti.ssm_daught_track1_score_mu_bck, "ssm_daught_track1_score_mu_bck/F");
    t_tagger->Branch("ssm_daught_track1_score_p_bck", &ti.ssm_daught_track1_score_p_bck, "ssm_daught_track1_score_p_bck/F");
    t_tagger->Branch("ssm_daught_track1_score_e_bck", &ti.ssm_daught_track1_score_e_bck, "ssm_daught_track1_score_e_bck/F");
    t_tagger->Branch("ssm_daught_track1_length", &ti.ssm_daught_track1_length, "ssm_daught_track1_length/F");
    t_tagger->Branch("ssm_daught_track1_direct_length", &ti.ssm_daught_track1_direct_length, "ssm_daught_track1_direct_length/F");
    t_tagger->Branch("ssm_daught_track1_length_ratio", &ti.ssm_daught_track1_length_ratio, "ssm_daught_track1_length_ratio/F");
    t_tagger->Branch("ssm_daught_track1_max_dev", &ti.ssm_daught_track1_max_dev, "ssm_daught_track1_max_dev/F");
    t_tagger->Branch("ssm_daught_track1_kine_energy_range", &ti.ssm_daught_track1_kine_energy_range, "ssm_daught_track1_kine_energy_range/F");
    t_tagger->Branch("ssm_daught_track1_kine_energy_range_mu", &ti.ssm_daught_track1_kine_energy_range_mu, "ssm_daught_track1_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_daught_track1_kine_energy_range_p", &ti.ssm_daught_track1_kine_energy_range_p, "ssm_daught_track1_kine_energy_range_p/F");
    t_tagger->Branch("ssm_daught_track1_kine_energy_range_e", &ti.ssm_daught_track1_kine_energy_range_e, "ssm_daught_track1_kine_energy_range_e/F");
    t_tagger->Branch("ssm_daught_track1_kine_energy_cal", &ti.ssm_daught_track1_kine_energy_cal, "ssm_daught_track1_kine_energy_cal/F");
    t_tagger->Branch("ssm_daught_track1_medium_dq_dx", &ti.ssm_daught_track1_medium_dq_dx, "ssm_daught_track1_medium_dq_dx/F");
    t_tagger->Branch("ssm_daught_track1_x_dir", &ti.ssm_daught_track1_x_dir, "ssm_daught_track1_x_dir/F");
    t_tagger->Branch("ssm_daught_track1_y_dir", &ti.ssm_daught_track1_y_dir, "ssm_daught_track1_y_dir/F");
    t_tagger->Branch("ssm_daught_track1_z_dir", &ti.ssm_daught_track1_z_dir, "ssm_daught_track1_z_dir/F");
    t_tagger->Branch("ssm_daught_track1_add_daught_track_counts_1", &ti.ssm_daught_track1_add_daught_track_counts_1, "ssm_daught_track1_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_daught_track1_add_daught_all_counts_1", &ti.ssm_daught_track1_add_daught_all_counts_1, "ssm_daught_track1_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_daught_track1_add_daught_track_counts_5", &ti.ssm_daught_track1_add_daught_track_counts_5, "ssm_daught_track1_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_daught_track1_add_daught_all_counts_5", &ti.ssm_daught_track1_add_daught_all_counts_5, "ssm_daught_track1_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_daught_track1_add_daught_track_counts_11", &ti.ssm_daught_track1_add_daught_track_counts_11, "ssm_daught_track1_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_daught_track1_add_daught_all_counts_11", &ti.ssm_daught_track1_add_daught_all_counts_11, "ssm_daught_track1_add_daught_all_counts_11/F");

    // SSM daughter track 2
    t_tagger->Branch("ssm_daught_track2_pdg", &ti.ssm_daught_track2_pdg, "ssm_daught_track2_pdg/F");
    t_tagger->Branch("ssm_daught_track2_score_mu_fwd", &ti.ssm_daught_track2_score_mu_fwd, "ssm_daught_track2_score_mu_fwd/F");
    t_tagger->Branch("ssm_daught_track2_score_p_fwd", &ti.ssm_daught_track2_score_p_fwd, "ssm_daught_track2_score_p_fwd/F");
    t_tagger->Branch("ssm_daught_track2_score_e_fwd", &ti.ssm_daught_track2_score_e_fwd, "ssm_daught_track2_score_e_fwd/F");
    t_tagger->Branch("ssm_daught_track2_score_mu_bck", &ti.ssm_daught_track2_score_mu_bck, "ssm_daught_track2_score_mu_bck/F");
    t_tagger->Branch("ssm_daught_track2_score_p_bck", &ti.ssm_daught_track2_score_p_bck, "ssm_daught_track2_score_p_bck/F");
    t_tagger->Branch("ssm_daught_track2_score_e_bck", &ti.ssm_daught_track2_score_e_bck, "ssm_daught_track2_score_e_bck/F");
    t_tagger->Branch("ssm_daught_track2_length", &ti.ssm_daught_track2_length, "ssm_daught_track2_length/F");
    t_tagger->Branch("ssm_daught_track2_direct_length", &ti.ssm_daught_track2_direct_length, "ssm_daught_track2_direct_length/F");
    t_tagger->Branch("ssm_daught_track2_length_ratio", &ti.ssm_daught_track2_length_ratio, "ssm_daught_track2_length_ratio/F");
    t_tagger->Branch("ssm_daught_track2_max_dev", &ti.ssm_daught_track2_max_dev, "ssm_daught_track2_max_dev/F");
    t_tagger->Branch("ssm_daught_track2_kine_energy_range", &ti.ssm_daught_track2_kine_energy_range, "ssm_daught_track2_kine_energy_range/F");
    t_tagger->Branch("ssm_daught_track2_kine_energy_range_mu", &ti.ssm_daught_track2_kine_energy_range_mu, "ssm_daught_track2_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_daught_track2_kine_energy_range_p", &ti.ssm_daught_track2_kine_energy_range_p, "ssm_daught_track2_kine_energy_range_p/F");
    t_tagger->Branch("ssm_daught_track2_kine_energy_range_e", &ti.ssm_daught_track2_kine_energy_range_e, "ssm_daught_track2_kine_energy_range_e/F");
    t_tagger->Branch("ssm_daught_track2_kine_energy_cal", &ti.ssm_daught_track2_kine_energy_cal, "ssm_daught_track2_kine_energy_cal/F");
    t_tagger->Branch("ssm_daught_track2_medium_dq_dx", &ti.ssm_daught_track2_medium_dq_dx, "ssm_daught_track2_medium_dq_dx/F");
    t_tagger->Branch("ssm_daught_track2_x_dir", &ti.ssm_daught_track2_x_dir, "ssm_daught_track2_x_dir/F");
    t_tagger->Branch("ssm_daught_track2_y_dir", &ti.ssm_daught_track2_y_dir, "ssm_daught_track2_y_dir/F");
    t_tagger->Branch("ssm_daught_track2_z_dir", &ti.ssm_daught_track2_z_dir, "ssm_daught_track2_z_dir/F");
    t_tagger->Branch("ssm_daught_track2_add_daught_track_counts_1", &ti.ssm_daught_track2_add_daught_track_counts_1, "ssm_daught_track2_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_daught_track2_add_daught_all_counts_1", &ti.ssm_daught_track2_add_daught_all_counts_1, "ssm_daught_track2_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_daught_track2_add_daught_track_counts_5", &ti.ssm_daught_track2_add_daught_track_counts_5, "ssm_daught_track2_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_daught_track2_add_daught_all_counts_5", &ti.ssm_daught_track2_add_daught_all_counts_5, "ssm_daught_track2_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_daught_track2_add_daught_track_counts_11", &ti.ssm_daught_track2_add_daught_track_counts_11, "ssm_daught_track2_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_daught_track2_add_daught_all_counts_11", &ti.ssm_daught_track2_add_daught_all_counts_11, "ssm_daught_track2_add_daught_all_counts_11/F");

    // SSM primary shower 1
    t_tagger->Branch("ssm_prim_shw1_pdg", &ti.ssm_prim_shw1_pdg, "ssm_prim_shw1_pdg/F");
    t_tagger->Branch("ssm_prim_shw1_score_mu_fwd", &ti.ssm_prim_shw1_score_mu_fwd, "ssm_prim_shw1_score_mu_fwd/F");
    t_tagger->Branch("ssm_prim_shw1_score_p_fwd", &ti.ssm_prim_shw1_score_p_fwd, "ssm_prim_shw1_score_p_fwd/F");
    t_tagger->Branch("ssm_prim_shw1_score_e_fwd", &ti.ssm_prim_shw1_score_e_fwd, "ssm_prim_shw1_score_e_fwd/F");
    t_tagger->Branch("ssm_prim_shw1_score_mu_bck", &ti.ssm_prim_shw1_score_mu_bck, "ssm_prim_shw1_score_mu_bck/F");
    t_tagger->Branch("ssm_prim_shw1_score_p_bck", &ti.ssm_prim_shw1_score_p_bck, "ssm_prim_shw1_score_p_bck/F");
    t_tagger->Branch("ssm_prim_shw1_score_e_bck", &ti.ssm_prim_shw1_score_e_bck, "ssm_prim_shw1_score_e_bck/F");
    t_tagger->Branch("ssm_prim_shw1_length", &ti.ssm_prim_shw1_length, "ssm_prim_shw1_length/F");
    t_tagger->Branch("ssm_prim_shw1_direct_length", &ti.ssm_prim_shw1_direct_length, "ssm_prim_shw1_direct_length/F");
    t_tagger->Branch("ssm_prim_shw1_length_ratio", &ti.ssm_prim_shw1_length_ratio, "ssm_prim_shw1_length_ratio/F");
    t_tagger->Branch("ssm_prim_shw1_max_dev", &ti.ssm_prim_shw1_max_dev, "ssm_prim_shw1_max_dev/F");
    t_tagger->Branch("ssm_prim_shw1_kine_energy_range", &ti.ssm_prim_shw1_kine_energy_range, "ssm_prim_shw1_kine_energy_range/F");
    t_tagger->Branch("ssm_prim_shw1_kine_energy_range_mu", &ti.ssm_prim_shw1_kine_energy_range_mu, "ssm_prim_shw1_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_prim_shw1_kine_energy_range_p", &ti.ssm_prim_shw1_kine_energy_range_p, "ssm_prim_shw1_kine_energy_range_p/F");
    t_tagger->Branch("ssm_prim_shw1_kine_energy_range_e", &ti.ssm_prim_shw1_kine_energy_range_e, "ssm_prim_shw1_kine_energy_range_e/F");
    t_tagger->Branch("ssm_prim_shw1_kine_energy_cal", &ti.ssm_prim_shw1_kine_energy_cal, "ssm_prim_shw1_kine_energy_cal/F");
    t_tagger->Branch("ssm_prim_shw1_kine_energy_best", &ti.ssm_prim_shw1_kine_energy_best, "ssm_prim_shw1_kine_energy_best/F");
    t_tagger->Branch("ssm_prim_shw1_medium_dq_dx", &ti.ssm_prim_shw1_medium_dq_dx, "ssm_prim_shw1_medium_dq_dx/F");
    t_tagger->Branch("ssm_prim_shw1_x_dir", &ti.ssm_prim_shw1_x_dir, "ssm_prim_shw1_x_dir/F");
    t_tagger->Branch("ssm_prim_shw1_y_dir", &ti.ssm_prim_shw1_y_dir, "ssm_prim_shw1_y_dir/F");
    t_tagger->Branch("ssm_prim_shw1_z_dir", &ti.ssm_prim_shw1_z_dir, "ssm_prim_shw1_z_dir/F");
    t_tagger->Branch("ssm_prim_shw1_add_daught_track_counts_1", &ti.ssm_prim_shw1_add_daught_track_counts_1, "ssm_prim_shw1_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_prim_shw1_add_daught_all_counts_1", &ti.ssm_prim_shw1_add_daught_all_counts_1, "ssm_prim_shw1_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_prim_shw1_add_daught_track_counts_5", &ti.ssm_prim_shw1_add_daught_track_counts_5, "ssm_prim_shw1_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_prim_shw1_add_daught_all_counts_5", &ti.ssm_prim_shw1_add_daught_all_counts_5, "ssm_prim_shw1_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_prim_shw1_add_daught_track_counts_11", &ti.ssm_prim_shw1_add_daught_track_counts_11, "ssm_prim_shw1_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_prim_shw1_add_daught_all_counts_11", &ti.ssm_prim_shw1_add_daught_all_counts_11, "ssm_prim_shw1_add_daught_all_counts_11/F");

    // SSM primary shower 2
    t_tagger->Branch("ssm_prim_shw2_pdg", &ti.ssm_prim_shw2_pdg, "ssm_prim_shw2_pdg/F");
    t_tagger->Branch("ssm_prim_shw2_score_mu_fwd", &ti.ssm_prim_shw2_score_mu_fwd, "ssm_prim_shw2_score_mu_fwd/F");
    t_tagger->Branch("ssm_prim_shw2_score_p_fwd", &ti.ssm_prim_shw2_score_p_fwd, "ssm_prim_shw2_score_p_fwd/F");
    t_tagger->Branch("ssm_prim_shw2_score_e_fwd", &ti.ssm_prim_shw2_score_e_fwd, "ssm_prim_shw2_score_e_fwd/F");
    t_tagger->Branch("ssm_prim_shw2_score_mu_bck", &ti.ssm_prim_shw2_score_mu_bck, "ssm_prim_shw2_score_mu_bck/F");
    t_tagger->Branch("ssm_prim_shw2_score_p_bck", &ti.ssm_prim_shw2_score_p_bck, "ssm_prim_shw2_score_p_bck/F");
    t_tagger->Branch("ssm_prim_shw2_score_e_bck", &ti.ssm_prim_shw2_score_e_bck, "ssm_prim_shw2_score_e_bck/F");
    t_tagger->Branch("ssm_prim_shw2_length", &ti.ssm_prim_shw2_length, "ssm_prim_shw2_length/F");
    t_tagger->Branch("ssm_prim_shw2_direct_length", &ti.ssm_prim_shw2_direct_length, "ssm_prim_shw2_direct_length/F");
    t_tagger->Branch("ssm_prim_shw2_length_ratio", &ti.ssm_prim_shw2_length_ratio, "ssm_prim_shw2_length_ratio/F");
    t_tagger->Branch("ssm_prim_shw2_max_dev", &ti.ssm_prim_shw2_max_dev, "ssm_prim_shw2_max_dev/F");
    t_tagger->Branch("ssm_prim_shw2_kine_energy_range", &ti.ssm_prim_shw2_kine_energy_range, "ssm_prim_shw2_kine_energy_range/F");
    t_tagger->Branch("ssm_prim_shw2_kine_energy_range_mu", &ti.ssm_prim_shw2_kine_energy_range_mu, "ssm_prim_shw2_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_prim_shw2_kine_energy_range_p", &ti.ssm_prim_shw2_kine_energy_range_p, "ssm_prim_shw2_kine_energy_range_p/F");
    t_tagger->Branch("ssm_prim_shw2_kine_energy_range_e", &ti.ssm_prim_shw2_kine_energy_range_e, "ssm_prim_shw2_kine_energy_range_e/F");
    t_tagger->Branch("ssm_prim_shw2_kine_energy_cal", &ti.ssm_prim_shw2_kine_energy_cal, "ssm_prim_shw2_kine_energy_cal/F");
    t_tagger->Branch("ssm_prim_shw2_kine_energy_best", &ti.ssm_prim_shw2_kine_energy_best, "ssm_prim_shw2_kine_energy_best/F");
    t_tagger->Branch("ssm_prim_shw2_medium_dq_dx", &ti.ssm_prim_shw2_medium_dq_dx, "ssm_prim_shw2_medium_dq_dx/F");
    t_tagger->Branch("ssm_prim_shw2_x_dir", &ti.ssm_prim_shw2_x_dir, "ssm_prim_shw2_x_dir/F");
    t_tagger->Branch("ssm_prim_shw2_y_dir", &ti.ssm_prim_shw2_y_dir, "ssm_prim_shw2_y_dir/F");
    t_tagger->Branch("ssm_prim_shw2_z_dir", &ti.ssm_prim_shw2_z_dir, "ssm_prim_shw2_z_dir/F");
    t_tagger->Branch("ssm_prim_shw2_add_daught_track_counts_1", &ti.ssm_prim_shw2_add_daught_track_counts_1, "ssm_prim_shw2_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_prim_shw2_add_daught_all_counts_1", &ti.ssm_prim_shw2_add_daught_all_counts_1, "ssm_prim_shw2_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_prim_shw2_add_daught_track_counts_5", &ti.ssm_prim_shw2_add_daught_track_counts_5, "ssm_prim_shw2_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_prim_shw2_add_daught_all_counts_5", &ti.ssm_prim_shw2_add_daught_all_counts_5, "ssm_prim_shw2_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_prim_shw2_add_daught_track_counts_11", &ti.ssm_prim_shw2_add_daught_track_counts_11, "ssm_prim_shw2_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_prim_shw2_add_daught_all_counts_11", &ti.ssm_prim_shw2_add_daught_all_counts_11, "ssm_prim_shw2_add_daught_all_counts_11/F");

    // SSM daughter shower 1
    t_tagger->Branch("ssm_daught_shw1_pdg", &ti.ssm_daught_shw1_pdg, "ssm_daught_shw1_pdg/F");
    t_tagger->Branch("ssm_daught_shw1_score_mu_fwd", &ti.ssm_daught_shw1_score_mu_fwd, "ssm_daught_shw1_score_mu_fwd/F");
    t_tagger->Branch("ssm_daught_shw1_score_p_fwd", &ti.ssm_daught_shw1_score_p_fwd, "ssm_daught_shw1_score_p_fwd/F");
    t_tagger->Branch("ssm_daught_shw1_score_e_fwd", &ti.ssm_daught_shw1_score_e_fwd, "ssm_daught_shw1_score_e_fwd/F");
    t_tagger->Branch("ssm_daught_shw1_score_mu_bck", &ti.ssm_daught_shw1_score_mu_bck, "ssm_daught_shw1_score_mu_bck/F");
    t_tagger->Branch("ssm_daught_shw1_score_p_bck", &ti.ssm_daught_shw1_score_p_bck, "ssm_daught_shw1_score_p_bck/F");
    t_tagger->Branch("ssm_daught_shw1_score_e_bck", &ti.ssm_daught_shw1_score_e_bck, "ssm_daught_shw1_score_e_bck/F");
    t_tagger->Branch("ssm_daught_shw1_length", &ti.ssm_daught_shw1_length, "ssm_daught_shw1_length/F");
    t_tagger->Branch("ssm_daught_shw1_direct_length", &ti.ssm_daught_shw1_direct_length, "ssm_daught_shw1_direct_length/F");
    t_tagger->Branch("ssm_daught_shw1_length_ratio", &ti.ssm_daught_shw1_length_ratio, "ssm_daught_shw1_length_ratio/F");
    t_tagger->Branch("ssm_daught_shw1_max_dev", &ti.ssm_daught_shw1_max_dev, "ssm_daught_shw1_max_dev/F");
    t_tagger->Branch("ssm_daught_shw1_kine_energy_range", &ti.ssm_daught_shw1_kine_energy_range, "ssm_daught_shw1_kine_energy_range/F");
    t_tagger->Branch("ssm_daught_shw1_kine_energy_range_mu", &ti.ssm_daught_shw1_kine_energy_range_mu, "ssm_daught_shw1_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_daught_shw1_kine_energy_range_p", &ti.ssm_daught_shw1_kine_energy_range_p, "ssm_daught_shw1_kine_energy_range_p/F");
    t_tagger->Branch("ssm_daught_shw1_kine_energy_range_e", &ti.ssm_daught_shw1_kine_energy_range_e, "ssm_daught_shw1_kine_energy_range_e/F");
    t_tagger->Branch("ssm_daught_shw1_kine_energy_cal", &ti.ssm_daught_shw1_kine_energy_cal, "ssm_daught_shw1_kine_energy_cal/F");
    t_tagger->Branch("ssm_daught_shw1_kine_energy_best", &ti.ssm_daught_shw1_kine_energy_best, "ssm_daught_shw1_kine_energy_best/F");
    t_tagger->Branch("ssm_daught_shw1_medium_dq_dx", &ti.ssm_daught_shw1_medium_dq_dx, "ssm_daught_shw1_medium_dq_dx/F");
    t_tagger->Branch("ssm_daught_shw1_x_dir", &ti.ssm_daught_shw1_x_dir, "ssm_daught_shw1_x_dir/F");
    t_tagger->Branch("ssm_daught_shw1_y_dir", &ti.ssm_daught_shw1_y_dir, "ssm_daught_shw1_y_dir/F");
    t_tagger->Branch("ssm_daught_shw1_z_dir", &ti.ssm_daught_shw1_z_dir, "ssm_daught_shw1_z_dir/F");
    t_tagger->Branch("ssm_daught_shw1_add_daught_track_counts_1", &ti.ssm_daught_shw1_add_daught_track_counts_1, "ssm_daught_shw1_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_daught_shw1_add_daught_all_counts_1", &ti.ssm_daught_shw1_add_daught_all_counts_1, "ssm_daught_shw1_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_daught_shw1_add_daught_track_counts_5", &ti.ssm_daught_shw1_add_daught_track_counts_5, "ssm_daught_shw1_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_daught_shw1_add_daught_all_counts_5", &ti.ssm_daught_shw1_add_daught_all_counts_5, "ssm_daught_shw1_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_daught_shw1_add_daught_track_counts_11", &ti.ssm_daught_shw1_add_daught_track_counts_11, "ssm_daught_shw1_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_daught_shw1_add_daught_all_counts_11", &ti.ssm_daught_shw1_add_daught_all_counts_11, "ssm_daught_shw1_add_daught_all_counts_11/F");

    // SSM daughter shower 2
    t_tagger->Branch("ssm_daught_shw2_pdg", &ti.ssm_daught_shw2_pdg, "ssm_daught_shw2_pdg/F");
    t_tagger->Branch("ssm_daught_shw2_score_mu_fwd", &ti.ssm_daught_shw2_score_mu_fwd, "ssm_daught_shw2_score_mu_fwd/F");
    t_tagger->Branch("ssm_daught_shw2_score_p_fwd", &ti.ssm_daught_shw2_score_p_fwd, "ssm_daught_shw2_score_p_fwd/F");
    t_tagger->Branch("ssm_daught_shw2_score_e_fwd", &ti.ssm_daught_shw2_score_e_fwd, "ssm_daught_shw2_score_e_fwd/F");
    t_tagger->Branch("ssm_daught_shw2_score_mu_bck", &ti.ssm_daught_shw2_score_mu_bck, "ssm_daught_shw2_score_mu_bck/F");
    t_tagger->Branch("ssm_daught_shw2_score_p_bck", &ti.ssm_daught_shw2_score_p_bck, "ssm_daught_shw2_score_p_bck/F");
    t_tagger->Branch("ssm_daught_shw2_score_e_bck", &ti.ssm_daught_shw2_score_e_bck, "ssm_daught_shw2_score_e_bck/F");
    t_tagger->Branch("ssm_daught_shw2_length", &ti.ssm_daught_shw2_length, "ssm_daught_shw2_length/F");
    t_tagger->Branch("ssm_daught_shw2_direct_length", &ti.ssm_daught_shw2_direct_length, "ssm_daught_shw2_direct_length/F");
    t_tagger->Branch("ssm_daught_shw2_length_ratio", &ti.ssm_daught_shw2_length_ratio, "ssm_daught_shw2_length_ratio/F");
    t_tagger->Branch("ssm_daught_shw2_max_dev", &ti.ssm_daught_shw2_max_dev, "ssm_daught_shw2_max_dev/F");
    t_tagger->Branch("ssm_daught_shw2_kine_energy_range", &ti.ssm_daught_shw2_kine_energy_range, "ssm_daught_shw2_kine_energy_range/F");
    t_tagger->Branch("ssm_daught_shw2_kine_energy_range_mu", &ti.ssm_daught_shw2_kine_energy_range_mu, "ssm_daught_shw2_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_daught_shw2_kine_energy_range_p", &ti.ssm_daught_shw2_kine_energy_range_p, "ssm_daught_shw2_kine_energy_range_p/F");
    t_tagger->Branch("ssm_daught_shw2_kine_energy_range_e", &ti.ssm_daught_shw2_kine_energy_range_e, "ssm_daught_shw2_kine_energy_range_e/F");
    t_tagger->Branch("ssm_daught_shw2_kine_energy_cal", &ti.ssm_daught_shw2_kine_energy_cal, "ssm_daught_shw2_kine_energy_cal/F");
    t_tagger->Branch("ssm_daught_shw2_kine_energy_best", &ti.ssm_daught_shw2_kine_energy_best, "ssm_daught_shw2_kine_energy_best/F");
    t_tagger->Branch("ssm_daught_shw2_medium_dq_dx", &ti.ssm_daught_shw2_medium_dq_dx, "ssm_daught_shw2_medium_dq_dx/F");
    t_tagger->Branch("ssm_daught_shw2_x_dir", &ti.ssm_daught_shw2_x_dir, "ssm_daught_shw2_x_dir/F");
    t_tagger->Branch("ssm_daught_shw2_y_dir", &ti.ssm_daught_shw2_y_dir, "ssm_daught_shw2_y_dir/F");
    t_tagger->Branch("ssm_daught_shw2_z_dir", &ti.ssm_daught_shw2_z_dir, "ssm_daught_shw2_z_dir/F");
    t_tagger->Branch("ssm_daught_shw2_add_daught_track_counts_1", &ti.ssm_daught_shw2_add_daught_track_counts_1, "ssm_daught_shw2_add_daught_track_counts_1/F");
    t_tagger->Branch("ssm_daught_shw2_add_daught_all_counts_1", &ti.ssm_daught_shw2_add_daught_all_counts_1, "ssm_daught_shw2_add_daught_all_counts_1/F");
    t_tagger->Branch("ssm_daught_shw2_add_daught_track_counts_5", &ti.ssm_daught_shw2_add_daught_track_counts_5, "ssm_daught_shw2_add_daught_track_counts_5/F");
    t_tagger->Branch("ssm_daught_shw2_add_daught_all_counts_5", &ti.ssm_daught_shw2_add_daught_all_counts_5, "ssm_daught_shw2_add_daught_all_counts_5/F");
    t_tagger->Branch("ssm_daught_shw2_add_daught_track_counts_11", &ti.ssm_daught_shw2_add_daught_track_counts_11, "ssm_daught_shw2_add_daught_track_counts_11/F");
    t_tagger->Branch("ssm_daught_shw2_add_daught_all_counts_11", &ti.ssm_daught_shw2_add_daught_all_counts_11, "ssm_daught_shw2_add_daught_all_counts_11/F");

    // SSM event-level angles and vertex
    t_tagger->Branch("ssm_nu_angle_z", &ti.ssm_nu_angle_z, "ssm_nu_angle_z/F");
    t_tagger->Branch("ssm_nu_angle_target", &ti.ssm_nu_angle_target, "ssm_nu_angle_target/F");
    t_tagger->Branch("ssm_nu_angle_absorber", &ti.ssm_nu_angle_absorber, "ssm_nu_angle_absorber/F");
    t_tagger->Branch("ssm_nu_angle_vertical", &ti.ssm_nu_angle_vertical, "ssm_nu_angle_vertical/F");
    t_tagger->Branch("ssm_con_nu_angle_z", &ti.ssm_con_nu_angle_z, "ssm_con_nu_angle_z/F");
    t_tagger->Branch("ssm_con_nu_angle_target", &ti.ssm_con_nu_angle_target, "ssm_con_nu_angle_target/F");
    t_tagger->Branch("ssm_con_nu_angle_absorber", &ti.ssm_con_nu_angle_absorber, "ssm_con_nu_angle_absorber/F");
    t_tagger->Branch("ssm_con_nu_angle_vertical", &ti.ssm_con_nu_angle_vertical, "ssm_con_nu_angle_vertical/F");
    t_tagger->Branch("ssm_prim_nu_angle_z", &ti.ssm_prim_nu_angle_z, "ssm_prim_nu_angle_z/F");
    t_tagger->Branch("ssm_prim_nu_angle_target", &ti.ssm_prim_nu_angle_target, "ssm_prim_nu_angle_target/F");
    t_tagger->Branch("ssm_prim_nu_angle_absorber", &ti.ssm_prim_nu_angle_absorber, "ssm_prim_nu_angle_absorber/F");
    t_tagger->Branch("ssm_prim_nu_angle_vertical", &ti.ssm_prim_nu_angle_vertical, "ssm_prim_nu_angle_vertical/F");
    t_tagger->Branch("ssm_track_angle_z", &ti.ssm_track_angle_z, "ssm_track_angle_z/F");
    t_tagger->Branch("ssm_track_angle_target", &ti.ssm_track_angle_target, "ssm_track_angle_target/F");
    t_tagger->Branch("ssm_track_angle_absorber", &ti.ssm_track_angle_absorber, "ssm_track_angle_absorber/F");
    t_tagger->Branch("ssm_track_angle_vertical", &ti.ssm_track_angle_vertical, "ssm_track_angle_vertical/F");
    t_tagger->Branch("ssm_vtxX", &ti.ssm_vtxX, "ssm_vtxX/F");
    t_tagger->Branch("ssm_vtxY", &ti.ssm_vtxY, "ssm_vtxY/F");
    t_tagger->Branch("ssm_vtxZ", &ti.ssm_vtxZ, "ssm_vtxZ/F");

    // SSM off-vertex activity
    t_tagger->Branch("ssm_offvtx_length", &ti.ssm_offvtx_length, "ssm_offvtx_length/F");
    t_tagger->Branch("ssm_offvtx_energy", &ti.ssm_offvtx_energy, "ssm_offvtx_energy/F");
    t_tagger->Branch("ssm_n_offvtx_tracks_1", &ti.ssm_n_offvtx_tracks_1, "ssm_n_offvtx_tracks_1/F");
    t_tagger->Branch("ssm_n_offvtx_tracks_3", &ti.ssm_n_offvtx_tracks_3, "ssm_n_offvtx_tracks_3/F");
    t_tagger->Branch("ssm_n_offvtx_tracks_5", &ti.ssm_n_offvtx_tracks_5, "ssm_n_offvtx_tracks_5/F");
    t_tagger->Branch("ssm_n_offvtx_tracks_8", &ti.ssm_n_offvtx_tracks_8, "ssm_n_offvtx_tracks_8/F");
    t_tagger->Branch("ssm_n_offvtx_tracks_11", &ti.ssm_n_offvtx_tracks_11, "ssm_n_offvtx_tracks_11/F");
    t_tagger->Branch("ssm_n_offvtx_showers_1", &ti.ssm_n_offvtx_showers_1, "ssm_n_offvtx_showers_1/F");
    t_tagger->Branch("ssm_n_offvtx_showers_3", &ti.ssm_n_offvtx_showers_3, "ssm_n_offvtx_showers_3/F");
    t_tagger->Branch("ssm_n_offvtx_showers_5", &ti.ssm_n_offvtx_showers_5, "ssm_n_offvtx_showers_5/F");
    t_tagger->Branch("ssm_n_offvtx_showers_8", &ti.ssm_n_offvtx_showers_8, "ssm_n_offvtx_showers_8/F");
    t_tagger->Branch("ssm_n_offvtx_showers_11", &ti.ssm_n_offvtx_showers_11, "ssm_n_offvtx_showers_11/F");

    // SSM off-vertex track 1
    t_tagger->Branch("ssm_offvtx_track1_pdg", &ti.ssm_offvtx_track1_pdg, "ssm_offvtx_track1_pdg/F");
    t_tagger->Branch("ssm_offvtx_track1_score_mu_fwd", &ti.ssm_offvtx_track1_score_mu_fwd, "ssm_offvtx_track1_score_mu_fwd/F");
    t_tagger->Branch("ssm_offvtx_track1_score_p_fwd", &ti.ssm_offvtx_track1_score_p_fwd, "ssm_offvtx_track1_score_p_fwd/F");
    t_tagger->Branch("ssm_offvtx_track1_score_e_fwd", &ti.ssm_offvtx_track1_score_e_fwd, "ssm_offvtx_track1_score_e_fwd/F");
    t_tagger->Branch("ssm_offvtx_track1_score_mu_bck", &ti.ssm_offvtx_track1_score_mu_bck, "ssm_offvtx_track1_score_mu_bck/F");
    t_tagger->Branch("ssm_offvtx_track1_score_p_bck", &ti.ssm_offvtx_track1_score_p_bck, "ssm_offvtx_track1_score_p_bck/F");
    t_tagger->Branch("ssm_offvtx_track1_score_e_bck", &ti.ssm_offvtx_track1_score_e_bck, "ssm_offvtx_track1_score_e_bck/F");
    t_tagger->Branch("ssm_offvtx_track1_length", &ti.ssm_offvtx_track1_length, "ssm_offvtx_track1_length/F");
    t_tagger->Branch("ssm_offvtx_track1_direct_length", &ti.ssm_offvtx_track1_direct_length, "ssm_offvtx_track1_direct_length/F");
    t_tagger->Branch("ssm_offvtx_track1_max_dev", &ti.ssm_offvtx_track1_max_dev, "ssm_offvtx_track1_max_dev/F");
    t_tagger->Branch("ssm_offvtx_track1_kine_energy_range", &ti.ssm_offvtx_track1_kine_energy_range, "ssm_offvtx_track1_kine_energy_range/F");
    t_tagger->Branch("ssm_offvtx_track1_kine_energy_range_mu", &ti.ssm_offvtx_track1_kine_energy_range_mu, "ssm_offvtx_track1_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_offvtx_track1_kine_energy_range_p", &ti.ssm_offvtx_track1_kine_energy_range_p, "ssm_offvtx_track1_kine_energy_range_p/F");
    t_tagger->Branch("ssm_offvtx_track1_kine_energy_range_e", &ti.ssm_offvtx_track1_kine_energy_range_e, "ssm_offvtx_track1_kine_energy_range_e/F");
    t_tagger->Branch("ssm_offvtx_track1_kine_energy_cal", &ti.ssm_offvtx_track1_kine_energy_cal, "ssm_offvtx_track1_kine_energy_cal/F");
    t_tagger->Branch("ssm_offvtx_track1_medium_dq_dx", &ti.ssm_offvtx_track1_medium_dq_dx, "ssm_offvtx_track1_medium_dq_dx/F");
    t_tagger->Branch("ssm_offvtx_track1_x_dir", &ti.ssm_offvtx_track1_x_dir, "ssm_offvtx_track1_x_dir/F");
    t_tagger->Branch("ssm_offvtx_track1_y_dir", &ti.ssm_offvtx_track1_y_dir, "ssm_offvtx_track1_y_dir/F");
    t_tagger->Branch("ssm_offvtx_track1_z_dir", &ti.ssm_offvtx_track1_z_dir, "ssm_offvtx_track1_z_dir/F");
    t_tagger->Branch("ssm_offvtx_track1_dist_mainvtx", &ti.ssm_offvtx_track1_dist_mainvtx, "ssm_offvtx_track1_dist_mainvtx/F");

    // SSM off-vertex shower 1
    t_tagger->Branch("ssm_offvtx_shw1_pdg_offvtx", &ti.ssm_offvtx_shw1_pdg_offvtx, "ssm_offvtx_shw1_pdg_offvtx/F");
    t_tagger->Branch("ssm_offvtx_shw1_score_mu_fwd", &ti.ssm_offvtx_shw1_score_mu_fwd, "ssm_offvtx_shw1_score_mu_fwd/F");
    t_tagger->Branch("ssm_offvtx_shw1_score_p_fwd", &ti.ssm_offvtx_shw1_score_p_fwd, "ssm_offvtx_shw1_score_p_fwd/F");
    t_tagger->Branch("ssm_offvtx_shw1_score_e_fwd", &ti.ssm_offvtx_shw1_score_e_fwd, "ssm_offvtx_shw1_score_e_fwd/F");
    t_tagger->Branch("ssm_offvtx_shw1_score_mu_bck", &ti.ssm_offvtx_shw1_score_mu_bck, "ssm_offvtx_shw1_score_mu_bck/F");
    t_tagger->Branch("ssm_offvtx_shw1_score_p_bck", &ti.ssm_offvtx_shw1_score_p_bck, "ssm_offvtx_shw1_score_p_bck/F");
    t_tagger->Branch("ssm_offvtx_shw1_score_e_bck", &ti.ssm_offvtx_shw1_score_e_bck, "ssm_offvtx_shw1_score_e_bck/F");
    t_tagger->Branch("ssm_offvtx_shw1_length", &ti.ssm_offvtx_shw1_length, "ssm_offvtx_shw1_length/F");
    t_tagger->Branch("ssm_offvtx_shw1_direct_length", &ti.ssm_offvtx_shw1_direct_length, "ssm_offvtx_shw1_direct_length/F");
    t_tagger->Branch("ssm_offvtx_shw1_max_dev", &ti.ssm_offvtx_shw1_max_dev, "ssm_offvtx_shw1_max_dev/F");
    t_tagger->Branch("ssm_offvtx_shw1_kine_energy_best", &ti.ssm_offvtx_shw1_kine_energy_best, "ssm_offvtx_shw1_kine_energy_best/F");
    t_tagger->Branch("ssm_offvtx_shw1_kine_energy_range", &ti.ssm_offvtx_shw1_kine_energy_range, "ssm_offvtx_shw1_kine_energy_range/F");
    t_tagger->Branch("ssm_offvtx_shw1_kine_energy_range_mu", &ti.ssm_offvtx_shw1_kine_energy_range_mu, "ssm_offvtx_shw1_kine_energy_range_mu/F");
    t_tagger->Branch("ssm_offvtx_shw1_kine_energy_range_p", &ti.ssm_offvtx_shw1_kine_energy_range_p, "ssm_offvtx_shw1_kine_energy_range_p/F");
    t_tagger->Branch("ssm_offvtx_shw1_kine_energy_range_e", &ti.ssm_offvtx_shw1_kine_energy_range_e, "ssm_offvtx_shw1_kine_energy_range_e/F");
    t_tagger->Branch("ssm_offvtx_shw1_kine_energy_cal", &ti.ssm_offvtx_shw1_kine_energy_cal, "ssm_offvtx_shw1_kine_energy_cal/F");
    t_tagger->Branch("ssm_offvtx_shw1_medium_dq_dx", &ti.ssm_offvtx_shw1_medium_dq_dx, "ssm_offvtx_shw1_medium_dq_dx/F");
    t_tagger->Branch("ssm_offvtx_shw1_x_dir", &ti.ssm_offvtx_shw1_x_dir, "ssm_offvtx_shw1_x_dir/F");
    t_tagger->Branch("ssm_offvtx_shw1_y_dir", &ti.ssm_offvtx_shw1_y_dir, "ssm_offvtx_shw1_y_dir/F");
    t_tagger->Branch("ssm_offvtx_shw1_z_dir", &ti.ssm_offvtx_shw1_z_dir, "ssm_offvtx_shw1_z_dir/F");
    t_tagger->Branch("ssm_offvtx_shw1_dist_mainvtx", &ti.ssm_offvtx_shw1_dist_mainvtx, "ssm_offvtx_shw1_dist_mainvtx/F");

    // SSM spacepoints (int scalars + vectors)
    t_tagger->Branch("ssmsp_Ntrack", &ti.ssmsp_Ntrack, "ssmsp_Ntrack/I");
    t_tagger->Branch("ssmsp_Nsp_tot", &ti.ssmsp_Nsp_tot, "ssmsp_Nsp_tot/I");
    t_tagger->Branch("ssmsp_Nsp", &ti.ssmsp_Nsp);
    t_tagger->Branch("ssmsp_pdg", &ti.ssmsp_pdg);
    t_tagger->Branch("ssmsp_id", &ti.ssmsp_id);
    t_tagger->Branch("ssmsp_mother", &ti.ssmsp_mother);
    t_tagger->Branch("ssmsp_x", &ti.ssmsp_x);
    t_tagger->Branch("ssmsp_y", &ti.ssmsp_y);
    t_tagger->Branch("ssmsp_z", &ti.ssmsp_z);
    t_tagger->Branch("ssmsp_dx", &ti.ssmsp_dx);
    t_tagger->Branch("ssmsp_dQ", &ti.ssmsp_dQ);
    t_tagger->Branch("ssmsp_KE", &ti.ssmsp_KE);
    t_tagger->Branch("ssmsp_containing_shower_id", &ti.ssmsp_containing_shower_id);
    t_tagger->Branch("ssmsp_containing_shower_ke", &ti.ssmsp_containing_shower_ke);
    t_tagger->Branch("ssmsp_containing_shower_flag", &ti.ssmsp_containing_shower_flag);

    // SSM kinematic variables
    t_tagger->Branch("ssm_kine_reco_Enu", &ti.ssm_kine_reco_Enu, "ssm_kine_reco_Enu/F");
    t_tagger->Branch("ssm_kine_reco_add_energy", &ti.ssm_kine_reco_add_energy, "ssm_kine_reco_add_energy/F");
    t_tagger->Branch("ssm_kine_energy_particle", &ti.ssm_kine_energy_particle);
    t_tagger->Branch("ssm_kine_energy_info", &ti.ssm_kine_energy_info);
    t_tagger->Branch("ssm_kine_particle_type", &ti.ssm_kine_particle_type);
    t_tagger->Branch("ssm_kine_energy_included", &ti.ssm_kine_energy_included);
    t_tagger->Branch("ssm_kine_pio_mass", &ti.ssm_kine_pio_mass, "ssm_kine_pio_mass/F");
    t_tagger->Branch("ssm_kine_pio_flag", &ti.ssm_kine_pio_flag, "ssm_kine_pio_flag/I");
    t_tagger->Branch("ssm_kine_pio_vtx_dis", &ti.ssm_kine_pio_vtx_dis, "ssm_kine_pio_vtx_dis/F");
    t_tagger->Branch("ssm_kine_pio_energy_1", &ti.ssm_kine_pio_energy_1, "ssm_kine_pio_energy_1/F");
    t_tagger->Branch("ssm_kine_pio_theta_1", &ti.ssm_kine_pio_theta_1, "ssm_kine_pio_theta_1/F");
    t_tagger->Branch("ssm_kine_pio_phi_1", &ti.ssm_kine_pio_phi_1, "ssm_kine_pio_phi_1/F");
    t_tagger->Branch("ssm_kine_pio_dis_1", &ti.ssm_kine_pio_dis_1, "ssm_kine_pio_dis_1/F");
    t_tagger->Branch("ssm_kine_pio_energy_2", &ti.ssm_kine_pio_energy_2, "ssm_kine_pio_energy_2/F");
    t_tagger->Branch("ssm_kine_pio_theta_2", &ti.ssm_kine_pio_theta_2, "ssm_kine_pio_theta_2/F");
    t_tagger->Branch("ssm_kine_pio_phi_2", &ti.ssm_kine_pio_phi_2, "ssm_kine_pio_phi_2/F");
    t_tagger->Branch("ssm_kine_pio_dis_2", &ti.ssm_kine_pio_dis_2, "ssm_kine_pio_dis_2/F");
    t_tagger->Branch("ssm_kine_pio_angle", &ti.ssm_kine_pio_angle, "ssm_kine_pio_angle/F");

    // ---- single-photon shower identification ----
    // (remaining TaggerInfo fields from shw_sp through photon_flag)
    // Due to the large number of branches (~700 remaining), they are registered
    // using the same pattern as above: field name == branch name, Float_t for
    // scalars, vector<float>/vector<int> for vectors.

    t_tagger->Branch("shw_sp_flag", &ti.shw_sp_flag, "shw_sp_flag/F");
    t_tagger->Branch("shw_sp_filled", &ti.shw_sp_filled, "shw_sp_filled/F");
    t_tagger->Branch("shw_sp_num_mip_tracks", &ti.shw_sp_num_mip_tracks, "shw_sp_num_mip_tracks/F");
    t_tagger->Branch("shw_sp_num_muons", &ti.shw_sp_num_muons, "shw_sp_num_muons/F");
    t_tagger->Branch("shw_sp_num_pions", &ti.shw_sp_num_pions, "shw_sp_num_pions/F");
    t_tagger->Branch("shw_sp_num_protons", &ti.shw_sp_num_protons, "shw_sp_num_protons/F");
    t_tagger->Branch("shw_sp_proton_length_1", &ti.shw_sp_proton_length_1, "shw_sp_proton_length_1/F");
    t_tagger->Branch("shw_sp_proton_dqdx_1", &ti.shw_sp_proton_dqdx_1, "shw_sp_proton_dqdx_1/F");
    t_tagger->Branch("shw_sp_proton_energy_1", &ti.shw_sp_proton_energy_1, "shw_sp_proton_energy_1/F");
    t_tagger->Branch("shw_sp_proton_length_2", &ti.shw_sp_proton_length_2, "shw_sp_proton_length_2/F");
    t_tagger->Branch("shw_sp_proton_dqdx_2", &ti.shw_sp_proton_dqdx_2, "shw_sp_proton_dqdx_2/F");
    t_tagger->Branch("shw_sp_proton_energy_2", &ti.shw_sp_proton_energy_2, "shw_sp_proton_energy_2/F");
    t_tagger->Branch("shw_sp_n_good_showers", &ti.shw_sp_n_good_showers, "shw_sp_n_good_showers/F");
    t_tagger->Branch("shw_sp_n_20mev_showers", &ti.shw_sp_n_20mev_showers, "shw_sp_n_20mev_showers/F");
    t_tagger->Branch("shw_sp_n_br1_showers", &ti.shw_sp_n_br1_showers, "shw_sp_n_br1_showers/F");
    t_tagger->Branch("shw_sp_n_br2_showers", &ti.shw_sp_n_br2_showers, "shw_sp_n_br2_showers/F");
    t_tagger->Branch("shw_sp_n_br3_showers", &ti.shw_sp_n_br3_showers, "shw_sp_n_br3_showers/F");
    t_tagger->Branch("shw_sp_n_br4_showers", &ti.shw_sp_n_br4_showers, "shw_sp_n_br4_showers/F");
    t_tagger->Branch("shw_sp_n_20br1_showers", &ti.shw_sp_n_20br1_showers, "shw_sp_n_20br1_showers/F");
    t_tagger->Branch("shw_sp_20mev_showers", &ti.shw_sp_20mev_showers);
    t_tagger->Branch("shw_sp_br1_showers", &ti.shw_sp_br1_showers);
    t_tagger->Branch("shw_sp_br2_showers", &ti.shw_sp_br2_showers);
    t_tagger->Branch("shw_sp_br3_showers", &ti.shw_sp_br3_showers);
    t_tagger->Branch("shw_sp_br4_showers", &ti.shw_sp_br4_showers);
    t_tagger->Branch("shw_sp_shw_vtx_dis", &ti.shw_sp_shw_vtx_dis, "shw_sp_shw_vtx_dis/F");
    t_tagger->Branch("shw_sp_max_shw_dis", &ti.shw_sp_max_shw_dis, "shw_sp_max_shw_dis/F");
    t_tagger->Branch("shw_sp_energy", &ti.shw_sp_energy, "shw_sp_energy/F");

    // shw_sp dQ/dx samples + MIP-like features
    t_tagger->Branch("shw_sp_vec_dQ_dx_0", &ti.shw_sp_vec_dQ_dx_0, "shw_sp_vec_dQ_dx_0/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_1", &ti.shw_sp_vec_dQ_dx_1, "shw_sp_vec_dQ_dx_1/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_2", &ti.shw_sp_vec_dQ_dx_2, "shw_sp_vec_dQ_dx_2/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_3", &ti.shw_sp_vec_dQ_dx_3, "shw_sp_vec_dQ_dx_3/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_4", &ti.shw_sp_vec_dQ_dx_4, "shw_sp_vec_dQ_dx_4/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_5", &ti.shw_sp_vec_dQ_dx_5, "shw_sp_vec_dQ_dx_5/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_6", &ti.shw_sp_vec_dQ_dx_6, "shw_sp_vec_dQ_dx_6/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_7", &ti.shw_sp_vec_dQ_dx_7, "shw_sp_vec_dQ_dx_7/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_8", &ti.shw_sp_vec_dQ_dx_8, "shw_sp_vec_dQ_dx_8/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_9", &ti.shw_sp_vec_dQ_dx_9, "shw_sp_vec_dQ_dx_9/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_10", &ti.shw_sp_vec_dQ_dx_10, "shw_sp_vec_dQ_dx_10/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_11", &ti.shw_sp_vec_dQ_dx_11, "shw_sp_vec_dQ_dx_11/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_12", &ti.shw_sp_vec_dQ_dx_12, "shw_sp_vec_dQ_dx_12/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_13", &ti.shw_sp_vec_dQ_dx_13, "shw_sp_vec_dQ_dx_13/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_14", &ti.shw_sp_vec_dQ_dx_14, "shw_sp_vec_dQ_dx_14/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_15", &ti.shw_sp_vec_dQ_dx_15, "shw_sp_vec_dQ_dx_15/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_16", &ti.shw_sp_vec_dQ_dx_16, "shw_sp_vec_dQ_dx_16/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_17", &ti.shw_sp_vec_dQ_dx_17, "shw_sp_vec_dQ_dx_17/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_18", &ti.shw_sp_vec_dQ_dx_18, "shw_sp_vec_dQ_dx_18/F");
    t_tagger->Branch("shw_sp_vec_dQ_dx_19", &ti.shw_sp_vec_dQ_dx_19, "shw_sp_vec_dQ_dx_19/F");
    t_tagger->Branch("shw_sp_max_dQ_dx_sample", &ti.shw_sp_max_dQ_dx_sample, "shw_sp_max_dQ_dx_sample/F");
    t_tagger->Branch("shw_sp_n_below_threshold", &ti.shw_sp_n_below_threshold, "shw_sp_n_below_threshold/F");
    t_tagger->Branch("shw_sp_n_below_zero", &ti.shw_sp_n_below_zero, "shw_sp_n_below_zero/F");
    t_tagger->Branch("shw_sp_n_lowest", &ti.shw_sp_n_lowest, "shw_sp_n_lowest/F");
    t_tagger->Branch("shw_sp_n_highest", &ti.shw_sp_n_highest, "shw_sp_n_highest/F");
    t_tagger->Branch("shw_sp_lowest_dQ_dx", &ti.shw_sp_lowest_dQ_dx, "shw_sp_lowest_dQ_dx/F");
    t_tagger->Branch("shw_sp_highest_dQ_dx", &ti.shw_sp_highest_dQ_dx, "shw_sp_highest_dQ_dx/F");
    t_tagger->Branch("shw_sp_medium_dQ_dx", &ti.shw_sp_medium_dQ_dx, "shw_sp_medium_dQ_dx/F");
    t_tagger->Branch("shw_sp_stem_length", &ti.shw_sp_stem_length, "shw_sp_stem_length/F");
    t_tagger->Branch("shw_sp_length_main", &ti.shw_sp_length_main, "shw_sp_length_main/F");
    t_tagger->Branch("shw_sp_length_total", &ti.shw_sp_length_total, "shw_sp_length_total/F");
    t_tagger->Branch("shw_sp_angle_beam", &ti.shw_sp_angle_beam, "shw_sp_angle_beam/F");
    t_tagger->Branch("shw_sp_iso_angle", &ti.shw_sp_iso_angle, "shw_sp_iso_angle/F");
    t_tagger->Branch("shw_sp_n_vertex", &ti.shw_sp_n_vertex, "shw_sp_n_vertex/F");
    t_tagger->Branch("shw_sp_n_good_tracks", &ti.shw_sp_n_good_tracks, "shw_sp_n_good_tracks/F");
    t_tagger->Branch("shw_sp_E_indirect_max_energy", &ti.shw_sp_E_indirect_max_energy, "shw_sp_E_indirect_max_energy/F");
    t_tagger->Branch("shw_sp_flag_all_above", &ti.shw_sp_flag_all_above, "shw_sp_flag_all_above/F");
    t_tagger->Branch("shw_sp_min_dQ_dx_5", &ti.shw_sp_min_dQ_dx_5, "shw_sp_min_dQ_dx_5/F");
    t_tagger->Branch("shw_sp_n_other_vertex", &ti.shw_sp_n_other_vertex, "shw_sp_n_other_vertex/F");
    t_tagger->Branch("shw_sp_n_stem_size", &ti.shw_sp_n_stem_size, "shw_sp_n_stem_size/F");
    t_tagger->Branch("shw_sp_flag_stem_trajectory", &ti.shw_sp_flag_stem_trajectory, "shw_sp_flag_stem_trajectory/F");
    t_tagger->Branch("shw_sp_min_dis", &ti.shw_sp_min_dis, "shw_sp_min_dis/F");
    t_tagger->Branch("shw_sp_vec_median_dedx", &ti.shw_sp_vec_median_dedx, "shw_sp_vec_median_dedx/F");
    t_tagger->Branch("shw_sp_vec_mean_dedx", &ti.shw_sp_vec_mean_dedx, "shw_sp_vec_mean_dedx/F");

    // shw_sp pi0
    t_tagger->Branch("shw_sp_pio_flag", &ti.shw_sp_pio_flag, "shw_sp_pio_flag/F");
    t_tagger->Branch("shw_sp_pio_mip_id", &ti.shw_sp_pio_mip_id, "shw_sp_pio_mip_id/F");
    t_tagger->Branch("shw_sp_pio_filled", &ti.shw_sp_pio_filled, "shw_sp_pio_filled/F");
    t_tagger->Branch("shw_sp_pio_flag_pio", &ti.shw_sp_pio_flag_pio, "shw_sp_pio_flag_pio/F");
    t_tagger->Branch("shw_sp_pio_1_flag", &ti.shw_sp_pio_1_flag, "shw_sp_pio_1_flag/F");
    t_tagger->Branch("shw_sp_pio_1_mass", &ti.shw_sp_pio_1_mass, "shw_sp_pio_1_mass/F");
    t_tagger->Branch("shw_sp_pio_1_pio_type", &ti.shw_sp_pio_1_pio_type, "shw_sp_pio_1_pio_type/F");
    t_tagger->Branch("shw_sp_pio_1_energy_1", &ti.shw_sp_pio_1_energy_1, "shw_sp_pio_1_energy_1/F");
    t_tagger->Branch("shw_sp_pio_1_energy_2", &ti.shw_sp_pio_1_energy_2, "shw_sp_pio_1_energy_2/F");
    t_tagger->Branch("shw_sp_pio_1_dis_1", &ti.shw_sp_pio_1_dis_1, "shw_sp_pio_1_dis_1/F");
    t_tagger->Branch("shw_sp_pio_1_dis_2", &ti.shw_sp_pio_1_dis_2, "shw_sp_pio_1_dis_2/F");
    t_tagger->Branch("shw_sp_pio_2_v_dis2", &ti.shw_sp_pio_2_v_dis2);
    t_tagger->Branch("shw_sp_pio_2_v_angle2", &ti.shw_sp_pio_2_v_angle2);
    t_tagger->Branch("shw_sp_pio_2_v_acc_length", &ti.shw_sp_pio_2_v_acc_length);
    t_tagger->Branch("shw_sp_pio_2_v_flag", &ti.shw_sp_pio_2_v_flag);

    // shw_sp low-energy michel
    t_tagger->Branch("shw_sp_lem_shower_total_length", &ti.shw_sp_lem_shower_total_length, "shw_sp_lem_shower_total_length/F");
    t_tagger->Branch("shw_sp_lem_shower_main_length", &ti.shw_sp_lem_shower_main_length, "shw_sp_lem_shower_main_length/F");
    t_tagger->Branch("shw_sp_lem_n_3seg", &ti.shw_sp_lem_n_3seg, "shw_sp_lem_n_3seg/F");
    t_tagger->Branch("shw_sp_lem_e_charge", &ti.shw_sp_lem_e_charge, "shw_sp_lem_e_charge/F");
    t_tagger->Branch("shw_sp_lem_e_dQdx", &ti.shw_sp_lem_e_dQdx, "shw_sp_lem_e_dQdx/F");
    t_tagger->Branch("shw_sp_lem_shower_num_segs", &ti.shw_sp_lem_shower_num_segs, "shw_sp_lem_shower_num_segs/F");
    t_tagger->Branch("shw_sp_lem_shower_num_main_segs", &ti.shw_sp_lem_shower_num_main_segs, "shw_sp_lem_shower_num_main_segs/F");
    t_tagger->Branch("shw_sp_lem_flag", &ti.shw_sp_lem_flag, "shw_sp_lem_flag/F");

    // shw_sp bad reconstruction br1-br4, lol, hol (all remaining shw_sp_* fields)
    // This is a very large block; I'll use a helper lambda to reduce visual noise.
    // Each line is: branch("name", &ti.name, "name/F") for scalars, branch("name", &ti.name) for vectors.

#define SCALAR_BR(name) t_tagger->Branch(#name, &ti.name, #name "/F")
#define VECTOR_BR(name) t_tagger->Branch(#name, &ti.name)

    SCALAR_BR(shw_sp_br_filled);
    SCALAR_BR(shw_sp_br1_flag); SCALAR_BR(shw_sp_br1_1_flag);
    SCALAR_BR(shw_sp_br1_1_shower_type); SCALAR_BR(shw_sp_br1_1_vtx_n_segs);
    SCALAR_BR(shw_sp_br1_1_energy); SCALAR_BR(shw_sp_br1_1_n_segs);
    SCALAR_BR(shw_sp_br1_1_flag_sg_topology); SCALAR_BR(shw_sp_br1_1_flag_sg_trajectory);
    SCALAR_BR(shw_sp_br1_1_sg_length);
    SCALAR_BR(shw_sp_br1_2_flag); SCALAR_BR(shw_sp_br1_2_energy);
    SCALAR_BR(shw_sp_br1_2_n_connected); SCALAR_BR(shw_sp_br1_2_max_length);
    SCALAR_BR(shw_sp_br1_2_n_connected_1); SCALAR_BR(shw_sp_br1_2_vtx_n_segs);
    SCALAR_BR(shw_sp_br1_2_n_shower_segs); SCALAR_BR(shw_sp_br1_2_max_length_ratio);
    SCALAR_BR(shw_sp_br1_2_shower_length);
    SCALAR_BR(shw_sp_br1_3_flag); SCALAR_BR(shw_sp_br1_3_energy);
    SCALAR_BR(shw_sp_br1_3_n_connected_p); SCALAR_BR(shw_sp_br1_3_max_length_p);
    SCALAR_BR(shw_sp_br1_3_n_shower_segs); SCALAR_BR(shw_sp_br1_3_flag_sg_topology);
    SCALAR_BR(shw_sp_br1_3_flag_sg_trajectory); SCALAR_BR(shw_sp_br1_3_n_shower_main_segs);
    SCALAR_BR(shw_sp_br1_3_sg_length);
    SCALAR_BR(shw_sp_br2_flag); SCALAR_BR(shw_sp_br2_flag_single_shower);
    SCALAR_BR(shw_sp_br2_num_valid_tracks); SCALAR_BR(shw_sp_br2_energy);
    SCALAR_BR(shw_sp_br2_angle1); SCALAR_BR(shw_sp_br2_angle2);
    SCALAR_BR(shw_sp_br2_angle); SCALAR_BR(shw_sp_br2_angle3);
    SCALAR_BR(shw_sp_br2_n_shower_main_segs); SCALAR_BR(shw_sp_br2_max_angle);
    SCALAR_BR(shw_sp_br2_sg_length); SCALAR_BR(shw_sp_br2_flag_sg_trajectory);
    SCALAR_BR(shw_sp_lol_flag); SCALAR_BR(shw_sp_lol_3_flag);
    SCALAR_BR(shw_sp_lol_3_angle_beam); SCALAR_BR(shw_sp_lol_3_min_angle);
    SCALAR_BR(shw_sp_lol_3_n_valid_tracks); SCALAR_BR(shw_sp_lol_3_vtx_n_segs);
    SCALAR_BR(shw_sp_lol_3_energy); SCALAR_BR(shw_sp_lol_3_shower_main_length);
    SCALAR_BR(shw_sp_lol_3_n_sum); SCALAR_BR(shw_sp_lol_3_n_out);
    VECTOR_BR(shw_sp_lol_1_v_energy); VECTOR_BR(shw_sp_lol_1_v_vtx_n_segs);
    VECTOR_BR(shw_sp_lol_1_v_nseg); VECTOR_BR(shw_sp_lol_1_v_angle);
    VECTOR_BR(shw_sp_lol_1_v_flag);
    VECTOR_BR(shw_sp_lol_2_v_length); VECTOR_BR(shw_sp_lol_2_v_angle);
    VECTOR_BR(shw_sp_lol_2_v_type); VECTOR_BR(shw_sp_lol_2_v_vtx_n_segs);
    VECTOR_BR(shw_sp_lol_2_v_energy); VECTOR_BR(shw_sp_lol_2_v_shower_main_length);
    VECTOR_BR(shw_sp_lol_2_v_flag_dir_weak); VECTOR_BR(shw_sp_lol_2_v_flag);

    // shw_sp br3
    SCALAR_BR(shw_sp_br3_1_energy); SCALAR_BR(shw_sp_br3_1_n_shower_segments);
    SCALAR_BR(shw_sp_br3_1_sg_flag_trajectory); SCALAR_BR(shw_sp_br3_1_sg_direct_length);
    SCALAR_BR(shw_sp_br3_1_sg_length); SCALAR_BR(shw_sp_br3_1_total_main_length);
    SCALAR_BR(shw_sp_br3_1_total_length); SCALAR_BR(shw_sp_br3_1_iso_angle);
    SCALAR_BR(shw_sp_br3_1_sg_flag_topology); SCALAR_BR(shw_sp_br3_1_flag);
    SCALAR_BR(shw_sp_br3_2_n_ele); SCALAR_BR(shw_sp_br3_2_n_other);
    SCALAR_BR(shw_sp_br3_2_energy); SCALAR_BR(shw_sp_br3_2_total_main_length);
    SCALAR_BR(shw_sp_br3_2_total_length); SCALAR_BR(shw_sp_br3_2_other_fid);
    SCALAR_BR(shw_sp_br3_2_flag);
    VECTOR_BR(shw_sp_br3_3_v_energy); VECTOR_BR(shw_sp_br3_3_v_angle);
    VECTOR_BR(shw_sp_br3_3_v_dir_length); VECTOR_BR(shw_sp_br3_3_v_length);
    VECTOR_BR(shw_sp_br3_3_v_flag);
    SCALAR_BR(shw_sp_br3_4_acc_length); SCALAR_BR(shw_sp_br3_4_total_length);
    SCALAR_BR(shw_sp_br3_4_energy); SCALAR_BR(shw_sp_br3_4_flag);
    VECTOR_BR(shw_sp_br3_5_v_dir_length); VECTOR_BR(shw_sp_br3_5_v_total_length);
    VECTOR_BR(shw_sp_br3_5_v_flag_avoid_muon_check); VECTOR_BR(shw_sp_br3_5_v_n_seg);
    VECTOR_BR(shw_sp_br3_5_v_angle); VECTOR_BR(shw_sp_br3_5_v_sg_length);
    VECTOR_BR(shw_sp_br3_5_v_energy); VECTOR_BR(shw_sp_br3_5_v_n_main_segs);
    VECTOR_BR(shw_sp_br3_5_v_n_segs); VECTOR_BR(shw_sp_br3_5_v_shower_main_length);
    VECTOR_BR(shw_sp_br3_5_v_shower_total_length); VECTOR_BR(shw_sp_br3_5_v_flag);
    VECTOR_BR(shw_sp_br3_6_v_angle); VECTOR_BR(shw_sp_br3_6_v_angle1);
    VECTOR_BR(shw_sp_br3_6_v_flag_shower_trajectory); VECTOR_BR(shw_sp_br3_6_v_direct_length);
    VECTOR_BR(shw_sp_br3_6_v_length); VECTOR_BR(shw_sp_br3_6_v_n_other_vtx_segs);
    VECTOR_BR(shw_sp_br3_6_v_energy); VECTOR_BR(shw_sp_br3_6_v_flag);
    SCALAR_BR(shw_sp_br3_7_energy); SCALAR_BR(shw_sp_br3_7_min_angle);
    SCALAR_BR(shw_sp_br3_7_sg_length); SCALAR_BR(shw_sp_br3_7_shower_main_length);
    SCALAR_BR(shw_sp_br3_7_flag);
    SCALAR_BR(shw_sp_br3_8_max_dQ_dx); SCALAR_BR(shw_sp_br3_8_energy);
    SCALAR_BR(shw_sp_br3_8_n_main_segs); SCALAR_BR(shw_sp_br3_8_shower_main_length);
    SCALAR_BR(shw_sp_br3_8_shower_length); SCALAR_BR(shw_sp_br3_8_flag);
    SCALAR_BR(shw_sp_br3_flag);

    // shw_sp br4
    SCALAR_BR(shw_sp_br4_1_shower_main_length); SCALAR_BR(shw_sp_br4_1_shower_total_length);
    SCALAR_BR(shw_sp_br4_1_min_dis); SCALAR_BR(shw_sp_br4_1_energy);
    SCALAR_BR(shw_sp_br4_1_flag_avoid_muon_check); SCALAR_BR(shw_sp_br4_1_n_vtx_segs);
    SCALAR_BR(shw_sp_br4_1_n_main_segs); SCALAR_BR(shw_sp_br4_1_flag);
    SCALAR_BR(shw_sp_br4_2_ratio_45); SCALAR_BR(shw_sp_br4_2_ratio_35);
    SCALAR_BR(shw_sp_br4_2_ratio_25); SCALAR_BR(shw_sp_br4_2_ratio_15);
    SCALAR_BR(shw_sp_br4_2_energy); SCALAR_BR(shw_sp_br4_2_ratio1_45);
    SCALAR_BR(shw_sp_br4_2_ratio1_35); SCALAR_BR(shw_sp_br4_2_ratio1_25);
    SCALAR_BR(shw_sp_br4_2_ratio1_15); SCALAR_BR(shw_sp_br4_2_iso_angle);
    SCALAR_BR(shw_sp_br4_2_iso_angle1); SCALAR_BR(shw_sp_br4_2_angle);
    SCALAR_BR(shw_sp_br4_2_flag); SCALAR_BR(shw_sp_br4_flag);

    // shw_sp hol
    SCALAR_BR(shw_sp_hol_1_n_valid_tracks); SCALAR_BR(shw_sp_hol_1_min_angle);
    SCALAR_BR(shw_sp_hol_1_energy); SCALAR_BR(shw_sp_hol_1_flag_all_shower);
    SCALAR_BR(shw_sp_hol_1_min_length); SCALAR_BR(shw_sp_hol_1_flag);
    SCALAR_BR(shw_sp_hol_2_min_angle); SCALAR_BR(shw_sp_hol_2_medium_dQ_dx);
    SCALAR_BR(shw_sp_hol_2_ncount); SCALAR_BR(shw_sp_hol_2_energy);
    SCALAR_BR(shw_sp_hol_2_flag); SCALAR_BR(shw_sp_hol_flag);

    // ---- pi0 identification ----
    SCALAR_BR(pio_flag); SCALAR_BR(pio_mip_id); SCALAR_BR(pio_filled);
    SCALAR_BR(pio_flag_pio); SCALAR_BR(pio_1_flag); SCALAR_BR(pio_1_mass);
    SCALAR_BR(pio_1_pio_type); SCALAR_BR(pio_1_energy_1); SCALAR_BR(pio_1_energy_2);
    SCALAR_BR(pio_1_dis_1); SCALAR_BR(pio_1_dis_2);
    VECTOR_BR(pio_2_v_dis2); VECTOR_BR(pio_2_v_angle2);
    VECTOR_BR(pio_2_v_acc_length); VECTOR_BR(pio_2_v_flag);

    // ---- single shower pi0 ----
    VECTOR_BR(sig_1_v_angle); VECTOR_BR(sig_1_v_flag_single_shower);
    VECTOR_BR(sig_1_v_energy); VECTOR_BR(sig_1_v_energy_1); VECTOR_BR(sig_1_v_flag);
    VECTOR_BR(sig_2_v_energy); VECTOR_BR(sig_2_v_shower_angle);
    VECTOR_BR(sig_2_v_flag_single_shower); VECTOR_BR(sig_2_v_medium_dQ_dx);
    VECTOR_BR(sig_2_v_start_dQ_dx); VECTOR_BR(sig_2_v_flag);
    SCALAR_BR(sig_flag);

    // ---- multiple gamma ----
    SCALAR_BR(mgo_energy); SCALAR_BR(mgo_max_energy); SCALAR_BR(mgo_total_energy);
    SCALAR_BR(mgo_n_showers); SCALAR_BR(mgo_max_energy_1); SCALAR_BR(mgo_max_energy_2);
    SCALAR_BR(mgo_total_other_energy); SCALAR_BR(mgo_n_total_showers);
    SCALAR_BR(mgo_total_other_energy_1); SCALAR_BR(mgo_flag);
    SCALAR_BR(mgt_flag_single_shower); SCALAR_BR(mgt_max_energy); SCALAR_BR(mgt_energy);
    SCALAR_BR(mgt_total_other_energy); SCALAR_BR(mgt_max_energy_1);
    SCALAR_BR(mgt_e_indirect_max_energy); SCALAR_BR(mgt_e_direct_max_energy);
    SCALAR_BR(mgt_n_direct_showers); SCALAR_BR(mgt_e_direct_total_energy);
    SCALAR_BR(mgt_flag_indirect_max_pio); SCALAR_BR(mgt_e_indirect_total_energy);
    SCALAR_BR(mgt_flag);

    // ---- shower to wall ----
    SCALAR_BR(stw_1_energy); SCALAR_BR(stw_1_dis); SCALAR_BR(stw_1_dQ_dx);
    SCALAR_BR(stw_1_flag_single_shower); SCALAR_BR(stw_1_n_pi0);
    SCALAR_BR(stw_1_num_valid_tracks); SCALAR_BR(stw_1_flag);
    VECTOR_BR(stw_2_v_medium_dQ_dx); VECTOR_BR(stw_2_v_energy);
    VECTOR_BR(stw_2_v_angle); VECTOR_BR(stw_2_v_dir_length);
    VECTOR_BR(stw_2_v_max_dQ_dx); VECTOR_BR(stw_2_v_flag);
    VECTOR_BR(stw_3_v_angle); VECTOR_BR(stw_3_v_dir_length);
    VECTOR_BR(stw_3_v_energy); VECTOR_BR(stw_3_v_medium_dQ_dx); VECTOR_BR(stw_3_v_flag);
    VECTOR_BR(stw_4_v_angle); VECTOR_BR(stw_4_v_dis);
    VECTOR_BR(stw_4_v_energy); VECTOR_BR(stw_4_v_flag);
    SCALAR_BR(stw_flag);

    // ---- single photon case ----
    SCALAR_BR(spt_flag_single_shower); SCALAR_BR(spt_energy);
    SCALAR_BR(spt_shower_main_length); SCALAR_BR(spt_shower_total_length);
    SCALAR_BR(spt_angle_beam); SCALAR_BR(spt_angle_vertical);
    SCALAR_BR(spt_max_dQ_dx); SCALAR_BR(spt_angle_beam_1);
    SCALAR_BR(spt_angle_drift); SCALAR_BR(spt_angle_drift_1);
    SCALAR_BR(spt_num_valid_tracks); SCALAR_BR(spt_n_vtx_segs);
    SCALAR_BR(spt_max_length); SCALAR_BR(spt_flag);

    // ---- stem length ----
    SCALAR_BR(stem_len_energy); SCALAR_BR(stem_len_length);
    SCALAR_BR(stem_len_flag_avoid_muon_check); SCALAR_BR(stem_len_num_daughters);
    SCALAR_BR(stem_len_daughter_length); SCALAR_BR(stem_len_flag);

    // ---- low-energy michel ----
    SCALAR_BR(lem_shower_total_length); SCALAR_BR(lem_shower_main_length);
    SCALAR_BR(lem_n_3seg); SCALAR_BR(lem_e_charge); SCALAR_BR(lem_e_dQdx);
    SCALAR_BR(lem_shower_num_segs); SCALAR_BR(lem_shower_num_main_segs); SCALAR_BR(lem_flag);

    // ---- broken muon ----
    SCALAR_BR(brm_n_mu_segs); SCALAR_BR(brm_Ep); SCALAR_BR(brm_energy);
    SCALAR_BR(brm_acc_length); SCALAR_BR(brm_shower_total_length);
    SCALAR_BR(brm_connected_length); SCALAR_BR(brm_n_size);
    SCALAR_BR(brm_acc_direct_length); SCALAR_BR(brm_n_shower_main_segs);
    SCALAR_BR(brm_n_mu_main); SCALAR_BR(brm_flag);

    // ---- compare muon energy ----
    SCALAR_BR(cme_mu_energy); SCALAR_BR(cme_energy); SCALAR_BR(cme_mu_length);
    SCALAR_BR(cme_length); SCALAR_BR(cme_angle_beam); SCALAR_BR(cme_flag);

    // ---- angular cut ----
    SCALAR_BR(anc_energy); SCALAR_BR(anc_angle); SCALAR_BR(anc_max_angle);
    SCALAR_BR(anc_max_length); SCALAR_BR(anc_acc_forward_length);
    SCALAR_BR(anc_acc_backward_length); SCALAR_BR(anc_acc_forward_length1);
    SCALAR_BR(anc_shower_main_length); SCALAR_BR(anc_shower_total_length);
    SCALAR_BR(anc_flag_main_outside); SCALAR_BR(anc_flag);

    // ---- stem direction ----
    SCALAR_BR(stem_dir_flag); SCALAR_BR(stem_dir_flag_single_shower);
    SCALAR_BR(stem_dir_filled); SCALAR_BR(stem_dir_angle); SCALAR_BR(stem_dir_energy);
    SCALAR_BR(stem_dir_angle1); SCALAR_BR(stem_dir_angle2); SCALAR_BR(stem_dir_angle3);
    SCALAR_BR(stem_dir_ratio);

    // ---- vertex inside shower ----
    SCALAR_BR(vis_1_filled); SCALAR_BR(vis_1_n_vtx_segs); SCALAR_BR(vis_1_energy);
    SCALAR_BR(vis_1_num_good_tracks); SCALAR_BR(vis_1_max_angle);
    SCALAR_BR(vis_1_max_shower_angle); SCALAR_BR(vis_1_tmp_length1);
    SCALAR_BR(vis_1_tmp_length2); SCALAR_BR(vis_1_particle_type); SCALAR_BR(vis_1_flag);
    SCALAR_BR(vis_2_filled); SCALAR_BR(vis_2_n_vtx_segs); SCALAR_BR(vis_2_min_angle);
    SCALAR_BR(vis_2_min_weak_track); SCALAR_BR(vis_2_angle_beam);
    SCALAR_BR(vis_2_min_angle1); SCALAR_BR(vis_2_iso_angle1);
    SCALAR_BR(vis_2_min_medium_dQ_dx); SCALAR_BR(vis_2_min_length);
    SCALAR_BR(vis_2_sg_length); SCALAR_BR(vis_2_max_angle);
    SCALAR_BR(vis_2_max_weak_track); SCALAR_BR(vis_2_flag);
    SCALAR_BR(vis_flag);

    // ---- bad reconstruction br1-br4 ----
    SCALAR_BR(br_filled); SCALAR_BR(br1_flag);
    SCALAR_BR(br1_1_flag); SCALAR_BR(br1_1_shower_type); SCALAR_BR(br1_1_vtx_n_segs);
    SCALAR_BR(br1_1_energy); SCALAR_BR(br1_1_n_segs);
    SCALAR_BR(br1_1_flag_sg_topology); SCALAR_BR(br1_1_flag_sg_trajectory);
    SCALAR_BR(br1_1_sg_length);
    SCALAR_BR(br1_2_flag); SCALAR_BR(br1_2_energy); SCALAR_BR(br1_2_n_connected);
    SCALAR_BR(br1_2_max_length); SCALAR_BR(br1_2_n_connected_1);
    SCALAR_BR(br1_2_vtx_n_segs); SCALAR_BR(br1_2_n_shower_segs);
    SCALAR_BR(br1_2_max_length_ratio); SCALAR_BR(br1_2_shower_length);
    SCALAR_BR(br1_3_flag); SCALAR_BR(br1_3_energy); SCALAR_BR(br1_3_n_connected_p);
    SCALAR_BR(br1_3_max_length_p); SCALAR_BR(br1_3_n_shower_segs);
    SCALAR_BR(br1_3_flag_sg_topology); SCALAR_BR(br1_3_flag_sg_trajectory);
    SCALAR_BR(br1_3_n_shower_main_segs); SCALAR_BR(br1_3_sg_length);

    SCALAR_BR(br2_flag); SCALAR_BR(br2_flag_single_shower); SCALAR_BR(br2_num_valid_tracks);
    SCALAR_BR(br2_energy); SCALAR_BR(br2_angle1); SCALAR_BR(br2_angle2);
    SCALAR_BR(br2_angle); SCALAR_BR(br2_angle3); SCALAR_BR(br2_n_shower_main_segs);
    SCALAR_BR(br2_max_angle); SCALAR_BR(br2_sg_length); SCALAR_BR(br2_flag_sg_trajectory);

    SCALAR_BR(br3_1_energy); SCALAR_BR(br3_1_n_shower_segments);
    SCALAR_BR(br3_1_sg_flag_trajectory); SCALAR_BR(br3_1_sg_direct_length);
    SCALAR_BR(br3_1_sg_length); SCALAR_BR(br3_1_total_main_length);
    SCALAR_BR(br3_1_total_length); SCALAR_BR(br3_1_iso_angle);
    SCALAR_BR(br3_1_sg_flag_topology); SCALAR_BR(br3_1_flag);
    SCALAR_BR(br3_2_n_ele); SCALAR_BR(br3_2_n_other); SCALAR_BR(br3_2_energy);
    SCALAR_BR(br3_2_total_main_length); SCALAR_BR(br3_2_total_length);
    SCALAR_BR(br3_2_other_fid); SCALAR_BR(br3_2_flag);
    VECTOR_BR(br3_3_v_energy); VECTOR_BR(br3_3_v_angle);
    VECTOR_BR(br3_3_v_dir_length); VECTOR_BR(br3_3_v_length); VECTOR_BR(br3_3_v_flag);
    SCALAR_BR(br3_4_acc_length); SCALAR_BR(br3_4_total_length);
    SCALAR_BR(br3_4_energy); SCALAR_BR(br3_4_flag);
    VECTOR_BR(br3_5_v_dir_length); VECTOR_BR(br3_5_v_total_length);
    VECTOR_BR(br3_5_v_flag_avoid_muon_check); VECTOR_BR(br3_5_v_n_seg);
    VECTOR_BR(br3_5_v_angle); VECTOR_BR(br3_5_v_sg_length);
    VECTOR_BR(br3_5_v_energy); VECTOR_BR(br3_5_v_n_main_segs);
    VECTOR_BR(br3_5_v_n_segs); VECTOR_BR(br3_5_v_shower_main_length);
    VECTOR_BR(br3_5_v_shower_total_length); VECTOR_BR(br3_5_v_flag);
    VECTOR_BR(br3_6_v_angle); VECTOR_BR(br3_6_v_angle1);
    VECTOR_BR(br3_6_v_flag_shower_trajectory); VECTOR_BR(br3_6_v_direct_length);
    VECTOR_BR(br3_6_v_length); VECTOR_BR(br3_6_v_n_other_vtx_segs);
    VECTOR_BR(br3_6_v_energy); VECTOR_BR(br3_6_v_flag);
    SCALAR_BR(br3_7_energy); SCALAR_BR(br3_7_min_angle); SCALAR_BR(br3_7_sg_length);
    SCALAR_BR(br3_7_shower_main_length); SCALAR_BR(br3_7_flag);
    SCALAR_BR(br3_8_max_dQ_dx); SCALAR_BR(br3_8_energy); SCALAR_BR(br3_8_n_main_segs);
    SCALAR_BR(br3_8_shower_main_length); SCALAR_BR(br3_8_shower_length); SCALAR_BR(br3_8_flag);
    SCALAR_BR(br3_flag);

    SCALAR_BR(br4_1_shower_main_length); SCALAR_BR(br4_1_shower_total_length);
    SCALAR_BR(br4_1_min_dis); SCALAR_BR(br4_1_energy);
    SCALAR_BR(br4_1_flag_avoid_muon_check); SCALAR_BR(br4_1_n_vtx_segs);
    SCALAR_BR(br4_1_n_main_segs); SCALAR_BR(br4_1_flag);
    SCALAR_BR(br4_2_ratio_45); SCALAR_BR(br4_2_ratio_35);
    SCALAR_BR(br4_2_ratio_25); SCALAR_BR(br4_2_ratio_15);
    SCALAR_BR(br4_2_energy); SCALAR_BR(br4_2_ratio1_45);
    SCALAR_BR(br4_2_ratio1_35); SCALAR_BR(br4_2_ratio1_25);
    SCALAR_BR(br4_2_ratio1_15); SCALAR_BR(br4_2_iso_angle);
    SCALAR_BR(br4_2_iso_angle1); SCALAR_BR(br4_2_angle);
    SCALAR_BR(br4_2_flag); SCALAR_BR(br4_flag);

    // ---- track overclustering ----
    VECTOR_BR(tro_1_v_particle_type); VECTOR_BR(tro_1_v_flag_dir_weak);
    VECTOR_BR(tro_1_v_min_dis); VECTOR_BR(tro_1_v_sg1_length);
    VECTOR_BR(tro_1_v_shower_main_length); VECTOR_BR(tro_1_v_max_n_vtx_segs);
    VECTOR_BR(tro_1_v_tmp_length); VECTOR_BR(tro_1_v_medium_dQ_dx);
    VECTOR_BR(tro_1_v_dQ_dx_cut); VECTOR_BR(tro_1_v_flag_shower_topology);
    VECTOR_BR(tro_1_v_flag);
    VECTOR_BR(tro_2_v_energy); VECTOR_BR(tro_2_v_stem_length);
    VECTOR_BR(tro_2_v_iso_angle); VECTOR_BR(tro_2_v_max_length);
    VECTOR_BR(tro_2_v_angle); VECTOR_BR(tro_2_v_flag);
    SCALAR_BR(tro_3_stem_length); SCALAR_BR(tro_3_n_muon_segs);
    SCALAR_BR(tro_3_energy); SCALAR_BR(tro_3_flag);
    VECTOR_BR(tro_4_v_dir2_mag); VECTOR_BR(tro_4_v_angle);
    VECTOR_BR(tro_4_v_angle1); VECTOR_BR(tro_4_v_angle2);
    VECTOR_BR(tro_4_v_length); VECTOR_BR(tro_4_v_length1);
    VECTOR_BR(tro_4_v_medium_dQ_dx); VECTOR_BR(tro_4_v_end_dQ_dx);
    VECTOR_BR(tro_4_v_energy); VECTOR_BR(tro_4_v_shower_main_length);
    VECTOR_BR(tro_4_v_flag_shower_trajectory); VECTOR_BR(tro_4_v_flag);
    VECTOR_BR(tro_5_v_max_angle); VECTOR_BR(tro_5_v_min_angle);
    VECTOR_BR(tro_5_v_max_length); VECTOR_BR(tro_5_v_iso_angle);
    VECTOR_BR(tro_5_v_n_vtx_segs); VECTOR_BR(tro_5_v_min_count);
    VECTOR_BR(tro_5_v_max_count); VECTOR_BR(tro_5_v_energy); VECTOR_BR(tro_5_v_flag);
    SCALAR_BR(tro_flag);

    // ---- high/low energy overlap ----
    SCALAR_BR(hol_1_n_valid_tracks); SCALAR_BR(hol_1_min_angle);
    SCALAR_BR(hol_1_energy); SCALAR_BR(hol_1_flag_all_shower);
    SCALAR_BR(hol_1_min_length); SCALAR_BR(hol_1_flag);
    SCALAR_BR(hol_2_min_angle); SCALAR_BR(hol_2_medium_dQ_dx);
    SCALAR_BR(hol_2_ncount); SCALAR_BR(hol_2_energy); SCALAR_BR(hol_2_flag);
    SCALAR_BR(hol_flag);
    SCALAR_BR(lol_flag);
    VECTOR_BR(lol_1_v_energy); VECTOR_BR(lol_1_v_vtx_n_segs);
    VECTOR_BR(lol_1_v_nseg); VECTOR_BR(lol_1_v_angle); VECTOR_BR(lol_1_v_flag);
    VECTOR_BR(lol_2_v_length); VECTOR_BR(lol_2_v_angle);
    VECTOR_BR(lol_2_v_type); VECTOR_BR(lol_2_v_vtx_n_segs);
    VECTOR_BR(lol_2_v_energy); VECTOR_BR(lol_2_v_shower_main_length);
    VECTOR_BR(lol_2_v_flag_dir_weak); VECTOR_BR(lol_2_v_flag);
    SCALAR_BR(lol_3_angle_beam); SCALAR_BR(lol_3_n_valid_tracks);
    SCALAR_BR(lol_3_min_angle); SCALAR_BR(lol_3_vtx_n_segs);
    SCALAR_BR(lol_3_energy); SCALAR_BR(lol_3_shower_main_length);
    SCALAR_BR(lol_3_n_out); SCALAR_BR(lol_3_n_sum); SCALAR_BR(lol_3_flag);

    // ---- cosmic tagger flags ----
    SCALAR_BR(cosmict_flag_1); SCALAR_BR(cosmict_flag_2); SCALAR_BR(cosmict_flag_3);
    SCALAR_BR(cosmict_flag_4); SCALAR_BR(cosmict_flag_5); SCALAR_BR(cosmict_flag_6);
    SCALAR_BR(cosmict_flag_7); SCALAR_BR(cosmict_flag_8); SCALAR_BR(cosmict_flag_9);
    VECTOR_BR(cosmict_flag_10);
    SCALAR_BR(cosmict_flag);

    // cosmict_2 (single muon)
    SCALAR_BR(cosmict_2_filled); SCALAR_BR(cosmict_2_particle_type);
    SCALAR_BR(cosmict_2_n_muon_tracks); SCALAR_BR(cosmict_2_total_shower_length);
    SCALAR_BR(cosmict_2_flag_inside); SCALAR_BR(cosmict_2_angle_beam);
    SCALAR_BR(cosmict_2_flag_dir_weak); SCALAR_BR(cosmict_2_dQ_dx_end);
    SCALAR_BR(cosmict_2_dQ_dx_front); SCALAR_BR(cosmict_2_theta);
    SCALAR_BR(cosmict_2_phi); SCALAR_BR(cosmict_2_valid_tracks);

    // cosmict_3 (single muon long)
    SCALAR_BR(cosmict_3_filled); SCALAR_BR(cosmict_3_flag_inside);
    SCALAR_BR(cosmict_3_angle_beam); SCALAR_BR(cosmict_3_flag_dir_weak);
    SCALAR_BR(cosmict_3_dQ_dx_end); SCALAR_BR(cosmict_3_dQ_dx_front);
    SCALAR_BR(cosmict_3_theta); SCALAR_BR(cosmict_3_phi); SCALAR_BR(cosmict_3_valid_tracks);

    // cosmict_4 (kinematics muon)
    SCALAR_BR(cosmict_4_filled); SCALAR_BR(cosmict_4_flag_inside);
    SCALAR_BR(cosmict_4_angle_beam); SCALAR_BR(cosmict_4_connected_showers);

    // cosmict_5 (kinematics muon long)
    SCALAR_BR(cosmict_5_filled); SCALAR_BR(cosmict_5_flag_inside);
    SCALAR_BR(cosmict_5_angle_beam); SCALAR_BR(cosmict_5_connected_showers);

    // cosmict_6 (special)
    SCALAR_BR(cosmict_6_filled); SCALAR_BR(cosmict_6_flag_dir_weak);
    SCALAR_BR(cosmict_6_flag_inside); SCALAR_BR(cosmict_6_angle);

    // cosmict_7 (muon + michel)
    SCALAR_BR(cosmict_7_filled); SCALAR_BR(cosmict_7_flag_sec);
    SCALAR_BR(cosmict_7_n_muon_tracks); SCALAR_BR(cosmict_7_total_shower_length);
    SCALAR_BR(cosmict_7_flag_inside); SCALAR_BR(cosmict_7_angle_beam);
    SCALAR_BR(cosmict_7_flag_dir_weak); SCALAR_BR(cosmict_7_dQ_dx_end);
    SCALAR_BR(cosmict_7_dQ_dx_front); SCALAR_BR(cosmict_7_theta); SCALAR_BR(cosmict_7_phi);

    // cosmict_8 (muon + michel + special)
    SCALAR_BR(cosmict_8_filled); SCALAR_BR(cosmict_8_flag_out);
    SCALAR_BR(cosmict_8_muon_length); SCALAR_BR(cosmict_8_acc_length);

    // cosmict_10 (front upstream dirt — vectors)
    VECTOR_BR(cosmict_10_flag_inside); VECTOR_BR(cosmict_10_vtx_z);
    VECTOR_BR(cosmict_10_flag_shower); VECTOR_BR(cosmict_10_flag_dir_weak);
    VECTOR_BR(cosmict_10_angle_beam); VECTOR_BR(cosmict_10_length);

    // ---- numu vs NC tagger ----
    SCALAR_BR(numu_cc_flag);
    VECTOR_BR(numu_cc_flag_1);
    VECTOR_BR(numu_cc_1_particle_type); VECTOR_BR(numu_cc_1_length);
    VECTOR_BR(numu_cc_1_medium_dQ_dx); VECTOR_BR(numu_cc_1_dQ_dx_cut);
    VECTOR_BR(numu_cc_1_direct_length); VECTOR_BR(numu_cc_1_n_daughter_tracks);
    VECTOR_BR(numu_cc_1_n_daughter_all);
    VECTOR_BR(numu_cc_flag_2);
    VECTOR_BR(numu_cc_2_length); VECTOR_BR(numu_cc_2_total_length);
    VECTOR_BR(numu_cc_2_n_daughter_tracks); VECTOR_BR(numu_cc_2_n_daughter_all);
    SCALAR_BR(numu_cc_flag_3);
    SCALAR_BR(numu_cc_3_particle_type); SCALAR_BR(numu_cc_3_max_length);
    SCALAR_BR(numu_cc_3_acc_track_length); SCALAR_BR(numu_cc_3_max_length_all);
    SCALAR_BR(numu_cc_3_max_muon_length); SCALAR_BR(numu_cc_3_n_daughter_tracks);
    SCALAR_BR(numu_cc_3_n_daughter_all);

    // ---- fiducial-volume flag ----
    SCALAR_BR(match_isFC);

    // ---- numu BDT scores ----
    SCALAR_BR(cosmict_2_4_score); SCALAR_BR(cosmict_3_5_score);
    SCALAR_BR(cosmict_6_score); SCALAR_BR(cosmict_7_score);
    SCALAR_BR(cosmict_8_score); SCALAR_BR(cosmict_10_score);
    SCALAR_BR(numu_1_score); SCALAR_BR(numu_2_score); SCALAR_BR(numu_3_score);
    SCALAR_BR(cosmict_score); SCALAR_BR(numu_score);

    // ---- nue BDT scores ----
    SCALAR_BR(mipid_score); SCALAR_BR(gap_score); SCALAR_BR(hol_lol_score);
    SCALAR_BR(cme_anc_score); SCALAR_BR(mgo_mgt_score); SCALAR_BR(br1_score);
    SCALAR_BR(br3_score); SCALAR_BR(br3_3_score); SCALAR_BR(br3_5_score);
    SCALAR_BR(br3_6_score); SCALAR_BR(stemdir_br2_score); SCALAR_BR(trimuon_score);
    SCALAR_BR(br4_tro_score); SCALAR_BR(mipquality_score);
    SCALAR_BR(pio_1_score); SCALAR_BR(pio_2_score); SCALAR_BR(stw_spt_score);
    SCALAR_BR(vis_1_score); SCALAR_BR(vis_2_score);
    SCALAR_BR(stw_2_score); SCALAR_BR(stw_3_score); SCALAR_BR(stw_4_score);
    SCALAR_BR(sig_1_score); SCALAR_BR(sig_2_score);
    SCALAR_BR(lol_1_score); SCALAR_BR(lol_2_score);
    SCALAR_BR(tro_1_score); SCALAR_BR(tro_2_score);
    SCALAR_BR(tro_4_score); SCALAR_BR(tro_5_score);
    SCALAR_BR(nue_score);

    SCALAR_BR(photon_flag);

#undef SCALAR_BR
#undef VECTOR_BR

    t_tagger->Fill();
    log->debug("UbooneTaggerOutputVisitor: wrote T_tagger with {} branches", t_tagger->GetNbranches());

    // ================================================================
    //  T_kine tree — all KineInfo fields
    // ================================================================
    TTree* t_kine = new TTree("T_kine", "T_kine");
    t_kine->SetDirectory(output_tf);

    t_kine->Branch("kine_nu_x_corr", &ki.kine_nu_x_corr, "kine_nu_x_corr/F");
    t_kine->Branch("kine_nu_y_corr", &ki.kine_nu_y_corr, "kine_nu_y_corr/F");
    t_kine->Branch("kine_nu_z_corr", &ki.kine_nu_z_corr, "kine_nu_z_corr/F");
    t_kine->Branch("kine_reco_Enu", &ki.kine_reco_Enu, "kine_reco_Enu/F");
    t_kine->Branch("kine_reco_add_energy", &ki.kine_reco_add_energy, "kine_reco_add_energy/F");
    t_kine->Branch("kine_energy_particle", &ki.kine_energy_particle);
    t_kine->Branch("kine_energy_info", &ki.kine_energy_info);
    t_kine->Branch("kine_particle_type", &ki.kine_particle_type);
    t_kine->Branch("kine_energy_included", &ki.kine_energy_included);
    t_kine->Branch("kine_pio_mass", &ki.kine_pio_mass, "kine_pio_mass/F");
    t_kine->Branch("kine_pio_flag", &ki.kine_pio_flag, "kine_pio_flag/I");
    t_kine->Branch("kine_pio_vtx_dis", &ki.kine_pio_vtx_dis, "kine_pio_vtx_dis/F");
    t_kine->Branch("kine_pio_energy_1", &ki.kine_pio_energy_1, "kine_pio_energy_1/F");
    t_kine->Branch("kine_pio_theta_1", &ki.kine_pio_theta_1, "kine_pio_theta_1/F");
    t_kine->Branch("kine_pio_phi_1", &ki.kine_pio_phi_1, "kine_pio_phi_1/F");
    t_kine->Branch("kine_pio_dis_1", &ki.kine_pio_dis_1, "kine_pio_dis_1/F");
    t_kine->Branch("kine_pio_energy_2", &ki.kine_pio_energy_2, "kine_pio_energy_2/F");
    t_kine->Branch("kine_pio_theta_2", &ki.kine_pio_theta_2, "kine_pio_theta_2/F");
    t_kine->Branch("kine_pio_phi_2", &ki.kine_pio_phi_2, "kine_pio_phi_2/F");
    t_kine->Branch("kine_pio_dis_2", &ki.kine_pio_dis_2, "kine_pio_dis_2/F");
    t_kine->Branch("kine_pio_angle", &ki.kine_pio_angle, "kine_pio_angle/F");

    t_kine->Fill();
    log->debug("UbooneTaggerOutputVisitor: wrote T_kine");

    output_tf->Write();
    output_tf->Close();
    delete output_tf;

    log->debug("UbooneTaggerOutputVisitor: updated {}", m_output_filename);
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
