#include "WireCellMatch/TimingTPCBundle.h"

#include <algorithm>
#include <cmath>
#include <iterator>

using namespace WireCell;
using namespace WireCell::Match;
using namespace WireCell::Clus::Facade;

namespace {
    // Kolmogorov-Smirnov statistic between two empirical CDFs constructed from
    // (already normalised) bin distributions.
    double calc_ks_test(const std::vector<double>& measured,
                        const std::vector<double>& predicted)
    {
        const std::size_t n = measured.size();
        std::vector<double> cum_m(n), cum_p(n);
        cum_m[0] = measured[0];
        cum_p[0] = predicted[0];
        for (std::size_t i = 1; i < n; ++i) {
            cum_m[i] = cum_m[i - 1] + measured[i];
            cum_p[i] = cum_p[i - 1] + predicted[i];
        }
        double max_diff = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            max_diff = std::max(max_diff, std::abs(cum_m[i] - cum_p[i]));
        }
        return max_diff;
    }
} // namespace

TimingTPCBundle::TimingTPCBundle(Opflash* flash, Cluster* main_cluster,
                                 int flash_index_id, int cluster_index_id)
    : flash(flash)
    , main_cluster(main_cluster)
    , orig_main_cluster(nullptr)
    , flash_index_id(flash_index_id)
    , cluster_index_id(cluster_index_id)
    , flag_close_to_PMT(false)
    , flag_at_x_boundary(false)
    , flag_spec_end(false)
    , flag_window_truncated(false)
    , flag_potential_bad_match(false)
    , flag_high_consistent(false)
    , ks_dis(1)
    , chi2(0)
    , ndf(0)
    , strength(0)
{
    // An unmatched cluster is represented by a bundle with a null flash (see
    // QLMatching: the "no matched flash" branch). Such a bundle carries no
    // predicted light, so there is no channel count to take from the flash.
    m_nchan = flash ? flash->get_num_channels() : 0;
    pred_flash.resize(m_nchan, 0);
    opdet_mask.resize(m_nchan, 0);
}

TimingTPCBundle::~TimingTPCBundle() = default;

double TimingTPCBundle::get_total_pred_light() const
{
    double sum = 0;
    for (auto v : pred_flash) sum += v;
    return sum;
}

bool TimingTPCBundle::examine_bundle(TimingTPCBundle* candidate_bundle)
{
    std::vector<double> measured_dist(m_nchan);
    std::vector<double> predicted_dist(m_nchan);
    std::vector<double> pe(m_nchan);
    std::vector<double> pe_err(m_nchan);
    std::vector<double> pred_pe(m_nchan);

    double candidate_ks_dis = 0;
    double candidate_chi2 = 0;

    for (int i = 0; i < m_nchan; ++i) {
        pe[i]     = flash->get_PE(i);
        pe_err[i] = flash->get_PE_err(i);
        pred_pe[i] = pred_flash.at(i) + candidate_bundle->get_pred_flash().at(i);
    }

    double total_predicted = 0;
    double total_measured  = 0;
    for (int j = 0; j < m_nchan; ++j) {
        measured_dist[j]  = pe[j];
        predicted_dist[j] = pred_pe[j];
        total_predicted  += pred_pe[j];
        total_measured   += pe[j];
    }
    if (total_predicted > 0) {
        for (int j = 0; j < m_nchan; ++j) predicted_dist[j] /= total_predicted;
    }
    if (total_measured > 0) {
        for (int j = 0; j < m_nchan; ++j) measured_dist[j] /= total_measured;
    }
    if (total_predicted > 0) candidate_ks_dis = calc_ks_test(measured_dist, predicted_dist);

    ndf = 0;
    for (int j = 0; j < m_nchan; ++j) {
        if (opdet_mask[j] == 0) continue;
        if (pe[j] < m_qp.pe_ndf_knee && pred_pe[j] < m_qp.pe_ndf_knee) { /* no-op */ }
        else ndf++;
        candidate_chi2 += std::pow(pred_pe[j] - pe[j], 2) / (pe[j] + std::pow(pe_err[j], 2));
    }

    if ((candidate_ks_dis < ks_dis || candidate_chi2 < chi2) &&
        (candidate_ks_dis < m_qp.ks_merge_max) && (candidate_chi2 / ndf < m_qp.chi2ndf_merge_max)) {
        return true;
    }
    return false;
}

void TimingTPCBundle::add_bundle(TimingTPCBundle* candidate_bundle)
{
    if (ks_dis * std::pow(chi2 / ndf, m_qp.addmerge_exponent) <
        candidate_bundle->get_ks_dis() *
            std::pow(candidate_bundle->get_chi2() / candidate_bundle->get_ndf(),
                     m_qp.addmerge_exponent)) {
        other_clusters.push_back(candidate_bundle->get_main_cluster());
        more_clusters.push_back(candidate_bundle->get_main_cluster());

        std::copy(candidate_bundle->get_other_clusters().begin(),
                  candidate_bundle->get_other_clusters().end(),
                  std::back_inserter(other_clusters));
        std::copy(candidate_bundle->get_more_clusters().begin(),
                  candidate_bundle->get_more_clusters().end(),
                  std::back_inserter(more_clusters));
    }
    else {
        other_clusters.push_back(main_cluster);
        more_clusters.push_back(main_cluster);

        main_cluster = candidate_bundle->get_main_cluster();
        std::copy(candidate_bundle->get_other_clusters().begin(),
                  candidate_bundle->get_other_clusters().end(),
                  std::back_inserter(other_clusters));
        std::copy(candidate_bundle->get_more_clusters().begin(),
                  candidate_bundle->get_more_clusters().end(),
                  std::back_inserter(more_clusters));

        flag_close_to_PMT  = candidate_bundle->get_flag_close_to_PMT();
        flag_at_x_boundary = candidate_bundle->get_flag_at_x_boundary();
        // Copy (not OR) to match the two flags above; revisit OR-vs-copy once a
        // consumer reads flag_window_truncated.
        flag_window_truncated = candidate_bundle->get_flag_window_truncated();
    }

    auto& pes = candidate_bundle->get_pred_flash();
    for (std::size_t i = 0; i < pred_flash.size(); ++i) pred_flash.at(i) += pes.at(i);
    examine_bundle();
}

bool TimingTPCBundle::examine_bundle()
{
    std::vector<double> measured_dist(m_nchan);
    std::vector<double> predicted_dist(m_nchan);
    std::vector<double> pe(m_nchan);
    std::vector<double> pe_err(m_nchan);
    std::vector<double> pred_pe(m_nchan);

    for (int i = 0; i < m_nchan; ++i) {
        pe[i]     = flash->get_PE(i);
        pe_err[i] = flash->get_PE_err(i);
        pred_pe[i] = pred_flash.at(i);
    }

    double total_predicted = 0;
    double total_measured  = 0;
    for (int j = 0; j < m_nchan; ++j) {
        measured_dist[j]  = pe[j];
        predicted_dist[j] = pred_pe[j];
        total_predicted  += pred_pe[j];
        total_measured   += pe[j];
    }
    if (total_predicted > 0) {
        for (int j = 0; j < m_nchan; ++j) predicted_dist[j] /= total_predicted;
    }
    if (total_measured > 0) {
        for (int j = 0; j < m_nchan; ++j) measured_dist[j] /= total_measured;
    }
    if (total_predicted > 0) ks_dis = calc_ks_test(measured_dist, predicted_dist);

    chi2 = 0;
    ndf  = 0;
    int nvalidopdets = 0;
    for (int j = 0; j < m_nchan; ++j) {
        if (opdet_mask[j] == 0) continue;
        ++nvalidopdets;
        if (pe[j] < m_qp.pe_ndf_knee && pred_pe[j] < m_qp.pe_ndf_knee) { /* no-op */ }
        else ndf++;
        chi2 += std::pow(pred_pe[j] - pe[j], 2) / (pe[j] + std::pow(pe_err[j], 2));
    }

    flag_high_consistent = false;
    if (ks_dis < m_qp.highconsist_ks_max && ndf >= m_qp.highconsist_min_ndf &&
        chi2 < ndf * nvalidopdets)
        flag_high_consistent = true;
    return flag_high_consistent;
}
