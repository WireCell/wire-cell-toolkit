#include "WireCellGen/BinnedDiffusion_transform.h"
#include "WireCellGen/GaussianDiffusion.h"
#include "WireCellUtil/Units.h"

#include <iostream>  // debug
#include <unordered_map>
using namespace std;

using namespace WireCell;

Gen::BinnedDiffusion_transform::BinnedDiffusion_transform(const Pimpos& pimpos, const Binning& tbins, double nsigma,
                                                          IRandom::pointer fluctuate,
                                                          ImpactDataCalculationStrategy calcstrat)
  : m_pimpos(pimpos)
  , m_tbins(tbins)
  , m_nsigma(nsigma)
  , m_fluctuate(fluctuate)
  , m_calcstrat(calcstrat)
  , m_window(0, 0)
  , m_outside_pitch(0)
  , m_outside_time(0)
{
}

bool Gen::BinnedDiffusion_transform::add(IDepo::pointer depo, double sigma_time, double sigma_pitch)
{
    const double center_time = depo->time();
    const double center_pitch = m_pimpos.distance(depo->pos());

    Gen::GausDesc time_desc(center_time, sigma_time);
    {
        double nmin_sigma = time_desc.distance(m_tbins.min());
        double nmax_sigma = time_desc.distance(m_tbins.max());

        double eff_nsigma = sigma_time > 0 ? m_nsigma : 0;
        if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) {
            // std::cerr << "BinnedDiffusion_transform: depo too far away in time sigma:"
            //           << " t_depo=" << center_time/units::ms << "ms not in:"
            //           << " t_bounds=[" << m_tbins.min()/units::ms << ","
            //           << m_tbins.max()/units::ms << "]ms"
            //           << " in Nsigma: [" << nmin_sigma << "," << nmax_sigma << "]\n";
            ++m_outside_time;
            return false;
        }
    }

    auto ibins = m_pimpos.impact_binning();

    Gen::GausDesc pitch_desc(center_pitch, sigma_pitch);
    {
        double nmin_sigma = pitch_desc.distance(ibins.min());
        double nmax_sigma = pitch_desc.distance(ibins.max());

        double eff_nsigma = sigma_pitch > 0 ? m_nsigma : 0;
        if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) {
            // std::cerr << "BinnedDiffusion_transform: depo too far away in pitch sigma: "
            //           << " p_depo=" << center_pitch/units::cm << "cm not in:"
            //           << " p_bounds=[" << ibins.min()/units::cm << ","
            //           << ibins.max()/units::cm << "]cm"
            //           << " in Nsigma:[" << nmin_sigma << "," << nmax_sigma << "]\n";
            ++m_outside_pitch;
            return false;
        }
    }

    // make GD and add to all covered impacts
    // int bin_beg = std::max(ibins.bin(center_pitch - sigma_pitch*m_nsigma), 0);
    // int bin_end = std::min(ibins.bin(center_pitch + sigma_pitch*m_nsigma)+1, ibins.nbins());
    // debug
    // int bin_center = ibins.bin(center_pitch);
    // cerr << "DEBUG center_pitch: "<<center_pitch/units::cm<<endl;
    // cerr << "DEBUG bin_center: "<<bin_center<<endl;

    auto gd = std::make_shared<GaussianDiffusion>(depo, time_desc, pitch_desc);
    // for (int bin = bin_beg; bin < bin_end; ++bin) {
    //   //   if (bin == bin_beg)  m_diffs.insert(gd);
    //   this->add(gd, bin);
    // }
    // m_diffs.insert(gd);
    m_diffs.push_back(gd);
    // std::cout << "yuhw m_diffs.size() " << m_diffs.size() << std::endl;
    return true;
}

// void Gen::BinnedDiffusion_transform::add(std::shared_ptr<GaussianDiffusion> gd, int bin)
// {
//     ImpactData::mutable_pointer idptr = nullptr;
//     auto it = m_impacts.find(bin);
//     if (it == m_impacts.end()) {
// 	idptr = std::make_shared<ImpactData>(bin);
// 	m_impacts[bin] = idptr;
//     }
//     else {
// 	idptr = it->second;
//     }
//     idptr->add(gd);
//     if (false) {                           // debug
//         auto mm = idptr->span();
//         cerr << "Gen::BinnedDiffusion_transform: add: "
//              << " poffoset="<<gd->poffset_bin()
//              << " toffoset="<<gd->toffset_bin()
//              << " charge=" << gd->depo()->charge()/units::eplus << " eles"
//              <<", for bin " << bin << " t=[" << mm.first/units::us << "," << mm.second/units::us << "]us\n";
//     }
//     m_diffs.insert(gd);
//     //m_diffs.push_back(gd);
// }

// void Gen::BinnedDiffusion_transform::erase(int begin_impact_number, int end_impact_number)
// {
//     for (int bin=begin_impact_number; bin<end_impact_number; ++bin) {
// 	m_impacts.erase(bin);
//     }
// }

void Gen::BinnedDiffusion_transform::get_charge_matrix(std::vector<Eigen::SparseMatrix<float>*>& vec_spmatrix,
                                                       std::vector<int>& vec_impact)
{
    const auto ib = m_pimpos.impact_binning();

    // map between reduced impact # to array #
    std::map<int, int> map_redimp_vec;
    for (size_t i = 0; i != vec_impact.size(); i++) {
        map_redimp_vec[vec_impact[i]] = int(i);
    }

    const auto rb = m_pimpos.region_binning();
    // map between impact # to channel #
    std::map<int, int> map_imp_ch;
    // map between impact # to reduced impact #
    std::map<int, int> map_imp_redimp;

    // std::cout << ib.nbins() << " " << rb.nbins() << std::endl;
    for (int wireind = 0; wireind != rb.nbins(); wireind++) {
        int wire_imp_no = m_pimpos.wire_impact(wireind);
        std::pair<int, int> imps_range = m_pimpos.wire_impacts(wireind);
        for (int imp_no = imps_range.first; imp_no != imps_range.second; imp_no++) {
            map_imp_ch[imp_no] = wireind;
            map_imp_redimp[imp_no] = imp_no - wire_imp_no;

            //  std::cout << imp_no << " " << wireind << " " << wire_imp_no << " " << ib.center(imp_no) << " " <<
            //  rb.center(wireind) << " " <<  ib.center(imp_no) - rb.center(wireind) << std::endl;
            // std::cout << imp_no << " " << map_imp_ch[imp_no] << " " << map_imp_redimp[imp_no] << std::endl;
        }
    }

    int min_imp = 0;
    int max_imp = ib.nbins();

    for (auto diff : m_diffs) {
        // std::cout << "yuhw m_diffs get_charge_matrix " << diff->depo()->time() << " " << diff->depo()->pos() << " " << diff->depo()->charge() / units::eplus
        //           << "e" << std::endl;
        // diff->set_sampling(m_tbins, ib, m_nsigma, 0, m_calcstrat);
        diff->set_sampling(m_tbins, ib, m_nsigma, m_fluctuate, m_calcstrat);
        // counter ++;

        const auto patch = diff->patch();
        const auto qweight = diff->weights();

        const int poffset_bin = diff->poffset_bin();
        const int toffset_bin = diff->toffset_bin();

        const int np = patch.rows();
        const int nt = patch.cols();

        for (int pbin = 0; pbin != np; pbin++) {
            int abs_pbin = pbin + poffset_bin;
            if (abs_pbin < min_imp || abs_pbin >= max_imp) continue;
            double weight = qweight[pbin];

            for (int tbin = 0; tbin != nt; tbin++) {
                int abs_tbin = tbin + toffset_bin;
                double charge = patch(pbin, tbin);

                // std::cout << map_redimp_vec[map_imp_redimp[abs_pbin] ] << " " <<
                // map_redimp_vec[map_imp_redimp[abs_pbin]+1] << " " << abs_tbin << " " << map_imp_ch[abs_pbin] <<
                // std::endl;

                vec_spmatrix.at(map_redimp_vec[map_imp_redimp[abs_pbin]])->coeffRef(abs_tbin, map_imp_ch[abs_pbin]) +=
                    charge * weight;
                vec_spmatrix.at(map_redimp_vec[map_imp_redimp[abs_pbin] + 1])
                    ->coeffRef(abs_tbin, map_imp_ch[abs_pbin]) += charge * (1 - weight);

                // if
                // (map_tuple_pos.find(std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]],map_imp_ch[abs_pbin],abs_tbin))==map_tuple_pos.end()){
                //   map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]],map_imp_ch[abs_pbin],abs_tbin)]
                //   = vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin] ]).size();
                //   vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]
                //   ]).push_back(std::make_tuple(map_imp_ch[abs_pbin],abs_tbin,charge*weight));
                // }else{
                //   std::get<2>(vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]
                //   ]).at(map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]],map_imp_ch[abs_pbin],abs_tbin)]))
                //   += charge * weight;
                // }

                // if
                // (map_tuple_pos.find(std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]+1],map_imp_ch[abs_pbin],abs_tbin))==map_tuple_pos.end()){
                //   map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]+1],map_imp_ch[abs_pbin],abs_tbin)]
                //   = vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]+1]).size();
                //   vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]+1]).push_back(std::make_tuple(map_imp_ch[abs_pbin],abs_tbin,charge*(1-weight)));
                // }else{
                //   std::get<2>(vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]+1]).at(map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]+1],map_imp_ch[abs_pbin],abs_tbin)])
                //   ) += charge*(1-weight);
                // }
            }
        }

        diff->clear_sampling();
        // need to figure out wire #, time #, charge, and weight ...
    }

    for (auto it = vec_spmatrix.begin(); it != vec_spmatrix.end(); it++) {
        (*it)->makeCompressed();
    }
}

// a new function to generate the result for the entire frame ...
void Gen::BinnedDiffusion_transform::get_charge_vec(
    std::vector<std::vector<std::tuple<int, int, double> > >& vec_vec_charge, std::vector<int>& vec_impact)
{
    const auto ib = m_pimpos.impact_binning();

    // map between reduced impact # to array #

    std::map<int, int> map_redimp_vec;
    std::vector<std::unordered_map<long int, int> > vec_map_pair_pos;
    for (size_t i = 0; i != vec_impact.size(); i++) {
        map_redimp_vec[vec_impact[i]] = int(i);
        std::unordered_map<long int, int> map_pair_pos;
        vec_map_pair_pos.push_back(map_pair_pos);
    }

    const auto rb = m_pimpos.region_binning();
    // map between impact # to channel #
    std::map<int, int> map_imp_ch;
    // map between impact # to reduced impact #
    std::map<int, int> map_imp_redimp;

    // std::cout << ib.nbins() << " " << rb.nbins() << std::endl;
    for (int wireind = 0; wireind != rb.nbins(); wireind++) {
        int wire_imp_no = m_pimpos.wire_impact(wireind);
        std::pair<int, int> imps_range = m_pimpos.wire_impacts(wireind);
        for (int imp_no = imps_range.first; imp_no != imps_range.second; imp_no++) {
            map_imp_ch[imp_no] = wireind;
            map_imp_redimp[imp_no] = imp_no - wire_imp_no;

            //  std::cout << imp_no << " " << wireind << " " << wire_imp_no << " " << ib.center(imp_no) << " " <<
            //  rb.center(wireind) << " " <<  ib.center(imp_no) - rb.center(wireind) << std::endl;
            // std::cout << imp_no << " " << map_imp_ch[imp_no] << " " << map_imp_redimp[imp_no] << std::endl;
        }
    }

    //  std::unordered_map<stgsd::tuple<int,int,int>, int> map_tuple_pos;

    //  int min_redimp =  m_pimpos.wire_impacts(2).first - m_pimpos.wire_impact(2);
    //  int max_redimp =  m_pimpos.wire_impacts(2).second - 1 - m_pimpos.wire_impact(2);
    int min_imp = 0;
    int max_imp = ib.nbins();

    // std::cout << min_redimp << " " << max_redimp << " " << max_imp << std::endl;

    int counter = 0;

    // std::set<std::shared_ptr<GaussianDiffusion>, GausDiffTimeCompare> m_diffs1;
    // for (auto diff : m_diffs){
    //    diff->set_sampling(m_tbins, ib, m_nsigma, m_fluctuate, m_calcstrat);
    //    m_diffs1.insert(diff);
    // }

    for (auto diff : m_diffs) {
        // std::cout << "yuhw m_diffs get_charge_vec " << diff->depo()->time() << " " << diff->depo()->pos() << " " << diff->depo()->charge() / units::eplus
        //           << "e" << std::endl;
        //    std::cout << diff->depo()->time() << std::endl
        // diff->set_sampling(m_tbins, ib, m_nsigma, 0, m_calcstrat);
        diff->set_sampling(m_tbins, ib, m_nsigma, m_fluctuate, m_calcstrat);
        counter++;

        const auto patch = diff->patch();
        const auto qweight = diff->weights();

        const int poffset_bin = diff->poffset_bin();
        const int toffset_bin = diff->toffset_bin();

        const int np = patch.rows();
        const int nt = patch.cols();

        // std::cout << np << " " << nt << std::endl;

        for (int pbin = 0; pbin != np; pbin++) {
            int abs_pbin = pbin + poffset_bin;
            if (abs_pbin < min_imp || abs_pbin >= max_imp) continue;
            double weight = qweight[pbin];
            auto const channel = map_imp_ch[abs_pbin];
            auto const redimp = map_imp_redimp[abs_pbin];
            auto const array_num_redimp = map_redimp_vec[redimp];
            auto const next_array_num_redimp = map_redimp_vec[redimp + 1];

            auto& map_pair_pos = vec_map_pair_pos.at(array_num_redimp);
            auto& next_map_pair_pos = vec_map_pair_pos.at(next_array_num_redimp);

            auto& vec_charge = vec_vec_charge.at(array_num_redimp);
            auto& next_vec_charge = vec_vec_charge.at(next_array_num_redimp);

            for (int tbin = 0; tbin != nt; tbin++) {
                int abs_tbin = tbin + toffset_bin;
                double charge = patch(pbin, tbin);

                // if (map_imp_ch[abs_pbin]==1459){
                //   std::cout << pbin+poffset_bin << " " << pbin << " " << tbin << " " << charge << " " << std::endl;
                // }
                //	std::cout << pbin << " " << tbin << " " << patch(pbin,tbin) << std::endl;
                // figure out how to convert the abs_pbin to fine position
                // figure out how to use the weight given the above ???
                // the other side
                //	if (map_imp_redimp[abs_pbin]==max_redimp){
                //  vec_vec_charge.at(map_redimp_vec[min_redimp]).push_back(std::make_tuple(map_imp_ch[abs_pbin]+1,abs_tbin,charge*(1-weight)));
                //}else{
                //}

                long int index1 = channel * 100000 + abs_tbin;
                auto it = map_pair_pos.find(index1);
                if (it == map_pair_pos.end()) {
                    map_pair_pos[index1] = vec_charge.size();
                    vec_charge.emplace_back(channel, abs_tbin, charge * weight);
                }
                else {
                    std::get<2>(vec_charge.at(it->second)) += charge * weight;
                }

                auto it1 = next_map_pair_pos.find(index1);
                if (it1 == next_map_pair_pos.end()) {
                    next_map_pair_pos[index1] = next_vec_charge.size();
                    next_vec_charge.emplace_back(channel, abs_tbin, charge * (1 - weight));
                }
                else {
                    std::get<2>(next_vec_charge.at(it1->second)) += charge * (1 - weight);
                }

                // if
                // (map_tuple_pos.find(std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]],map_imp_ch[abs_pbin],abs_tbin))==map_tuple_pos.end()){
                //   map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]],map_imp_ch[abs_pbin],abs_tbin)]
                //   = vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin] ]).size();
                //   vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]
                //   ]).push_back(std::make_tuple(map_imp_ch[abs_pbin],abs_tbin,charge*weight));
                // }else{
                //   std::get<2>(vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]
                //   ]).at(map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]],map_imp_ch[abs_pbin],abs_tbin)]))
                //   += charge * weight;
                // }

                // if
                // (map_tuple_pos.find(std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]+1],map_imp_ch[abs_pbin],abs_tbin))==map_tuple_pos.end()){
                //   map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]+1],map_imp_ch[abs_pbin],abs_tbin)]
                //   = vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]+1]).size();
                //   vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]+1]).push_back(std::make_tuple(map_imp_ch[abs_pbin],abs_tbin,charge*(1-weight)));
                // }else{
                //   std::get<2>(vec_vec_charge.at(map_redimp_vec[map_imp_redimp[abs_pbin]+1]).at(map_tuple_pos[std::make_tuple(map_redimp_vec[map_imp_redimp[abs_pbin]+1],map_imp_ch[abs_pbin],abs_tbin)])
                //   ) += charge*(1-weight);
                // }
            }
        }

        if (counter % 5000 == 0) {
            // std::vector<std::tuple<int,int,int> > del_keys;
            // for (auto it  = map_tuple_pos.begin(); it!=map_tuple_pos.end(); it++){
            // 	if (get<2>(it->first) < toffset_bin - 60){
            // 	  del_keys.push_back(it->first);
            // 	}
            // }
            // for (auto it = del_keys.begin(); it!=del_keys.end(); it++){
            // 	map_tuple_pos.erase(*it);
            // }
            // map_tuple_pos.clear();
            for (auto it = vec_map_pair_pos.begin(); it != vec_map_pair_pos.end(); it++) {
                it->clear();
            }
        }

        diff->clear_sampling();
        // need to figure out wire #, time #, charge, and weight ...
    }

    //
}

// Gen::ImpactData::pointer Gen::BinnedDiffusion_transform::impact_data(int bin) const
// {
//     const auto ib = m_pimpos.impact_binning();
//     if (! ib.inbounds(bin)) {
//         return nullptr;
//     }

//     auto it = m_impacts.find(bin);
//     if (it == m_impacts.end()) {
// 	return nullptr;
//     }
//     auto idptr = it->second;

//     // make sure all diffusions have been sampled
//     for (auto diff : idptr->diffusions()) {
//       diff->set_sampling(m_tbins, ib, m_nsigma, m_fluctuate, m_calcstrat);
//       //diff->set_sampling(m_tbins, ib, m_nsigma, 0, m_calcstrat);
//     }

//     idptr->calculate(m_tbins.nbins());
//     return idptr;
// }

static std::pair<double, double> gausdesc_range(const std::vector<Gen::GausDesc> gds, double nsigma)
{
    int ncount = -1;
    double vmin = 0, vmax = 0;
    for (auto gd : gds) {
        ++ncount;

        const double lvmin = gd.center - gd.sigma * nsigma;
        const double lvmax = gd.center + gd.sigma * nsigma;
        if (!ncount) {
            vmin = lvmin;
            vmax = lvmax;
            continue;
        }
        vmin = std::min(vmin, lvmin);
        vmax = std::max(vmax, lvmax);
    }
    return std::make_pair(vmin, vmax);
}

std::pair<double, double> Gen::BinnedDiffusion_transform::pitch_range(double nsigma) const
{
    std::vector<Gen::GausDesc> gds;
    for (auto diff : m_diffs) {
        gds.push_back(diff->pitch_desc());
    }
    return gausdesc_range(gds, nsigma);
}

std::pair<int, int> Gen::BinnedDiffusion_transform::impact_bin_range(double nsigma) const
{
    const auto ibins = m_pimpos.impact_binning();
    auto mm = pitch_range(nsigma);
    return std::make_pair(std::max(ibins.bin(mm.first), 0), std::min(ibins.bin(mm.second) + 1, ibins.nbins()));
}

std::pair<double, double> Gen::BinnedDiffusion_transform::time_range(double nsigma) const
{
    std::vector<Gen::GausDesc> gds;
    for (auto diff : m_diffs) {
        gds.push_back(diff->time_desc());
    }
    return gausdesc_range(gds, nsigma);
}

std::pair<int, int> Gen::BinnedDiffusion_transform::time_bin_range(double nsigma) const
{
    auto mm = time_range(nsigma);
    return std::make_pair(std::max(m_tbins.bin(mm.first), 0), std::min(m_tbins.bin(mm.second) + 1, m_tbins.nbins()));
}
