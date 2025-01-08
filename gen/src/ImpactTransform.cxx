#include "WireCellGen/ImpactTransform.h"

#include "WireCellAux/DftTools.h"

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellUtil/Exceptions.h"

#include <iostream>  // debugging.
using namespace std;

using namespace WireCell;
using WireCell::Aux::DftTools::fwd;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv;
using WireCell::Aux::DftTools::inv_c2r;


Gen::ImpactTransform::ImpactTransform(IPlaneImpactResponse::pointer pir,
                                      const IDFT::pointer& dft,
                                      BinnedDiffusion_transform& bd)
  : m_pir(pir)
  , m_dft(dft)
  , m_bd(bd)
{
    // arrange the field response (210 in total, pitch_range/impact)
    // number of wires nwires ...
    m_num_group = std::round(m_pir->pitch() / m_pir->impact()) + 1;  // 11
    m_num_pad_wire = std::round((m_pir->nwires() - 1) / 2.);         // 10 for wires, 5 for PCB strips

    // const auto pimpos = m_bd.pimpos();
    
    //std::cerr << "ImpactTransform: num_group:" << m_num_group << " num_pad_wire:" << m_num_pad_wire << std::endl;
    for (int i = 0; i != m_num_group; i++) {
        double rel_cen_imp_pos;
        if (i != m_num_group - 1) {
            rel_cen_imp_pos = -m_pir->pitch() / 2. + m_pir->impact() * i + 1e-9;
        }
        else {
            rel_cen_imp_pos = -m_pir->pitch() / 2. + m_pir->impact() * i - 1e-9;
        }
        m_vec_impact.push_back(std::round(rel_cen_imp_pos / m_pir->impact()));
        std::map<int, IImpactResponse::pointer> map_resp;  // already in freq domain

        //std::cerr << "ImpactTransform: " << rel_cen_imp_pos << std::endl;
        for (int j = 0; j != m_pir->nwires(); j++) {

            try {
                map_resp[j - m_num_pad_wire] = m_pir->closest(rel_cen_imp_pos - (j - m_num_pad_wire) * m_pir->pitch());
            }
            catch (ValueError& ve) {
                std::cerr << "ImpactTransform: I angered PIR with: i="
                          << i << " j=" << j
                          << " rel_cen_imp_pos=" << rel_cen_imp_pos
                          << " m_num_pad_wire=" << m_num_pad_wire
                          << " m_num_group=" << m_num_group
                          << " pir->pitch=" << m_pir->pitch()
                          << " pir->impact=" << m_pir->impact()
                          << " pir->nwires=" << m_pir->nwires()
                          << " looking for: "
                          << rel_cen_imp_pos - (j - m_num_pad_wire) * m_pir->pitch()
                          << std::endl;                
                continue;
            }
                
            Waveform::compseq_t response_spectrum = map_resp[j - m_num_pad_wire]->spectrum();

        }

        m_vec_map_resp.push_back(map_resp);

        std::vector<std::tuple<int, int, double> > vec_charge;  // ch, time, charge
        m_vec_vec_charge.push_back(vec_charge);
    }

    // now work on the charge part ...
    // trying to sampling ...
    m_bd.get_charge_vec(m_vec_vec_charge, m_vec_impact);
    // std::cout << nwires << " " << nsamples << std::endl;

    // length and width ...

    //    std::cout << nwires << " " << nsamples << std::endl;
    std::pair<int, int> impact_range = m_bd.impact_bin_range(m_bd.get_nsigma());
    std::pair<int, int> time_range = m_bd.time_bin_range(m_bd.get_nsigma());

    //  std::cout << impact_range.first << " " << impact_range.second << " " << time_range.first << " " <<
    //  time_range.second << std::endl;

    int start_ch = std::floor(impact_range.first * 1.0 / (m_num_group - 1)) - 1;
    int end_ch = std::ceil(impact_range.second * 1.0 / (m_num_group - 1)) + 2;
    if ((end_ch - start_ch) % 2 == 1) end_ch += 1;
    int start_tick = time_range.first - 1;
    int end_tick = time_range.second + 2;
    if ((end_tick - start_tick) % 2 == 1) end_tick += 1;

    //  m_decon_data = Array::array_xxf::Zero(nwires+2*m_num_pad_wire,nsamples);
    // for saving the accumulated wire data in the time frequency domain ...
    // adding no padding now, it make the FFT slower, need some other methods ...

    int npad_wire = 0;
    const size_t ntotal_wires = fft_best_length(end_ch - start_ch + 2 * m_num_pad_wire, 1);

    npad_wire = (ntotal_wires - end_ch + start_ch) / 2;
    m_start_ch = start_ch - npad_wire;
    m_end_ch = end_ch + npad_wire;
    // std::cout << start_ch << " " << end_ch << " " << npad_wire << " " << start_tick << " " << end_tick << " " <<
    // m_start_ch << " " << m_end_ch << std::endl;

    int npad_time = m_pir->closest(0)->waveform_pad();
    const size_t ntotal_ticks = fft_best_length(end_tick - start_tick + npad_time);

    npad_time = ntotal_ticks - end_tick + start_tick;
    m_start_tick = start_tick;
    m_end_tick = end_tick + npad_time;

    Array::array_xxc acc_data_f_w =
        Array::array_xxc::Zero(end_ch - start_ch + 2 * npad_wire, m_end_tick - m_start_tick);

    int num_double = (m_vec_vec_charge.size() - 1) / 2;

    // speed up version , first five
    for (int i = 0; i != num_double; i++) {
        // if (i!=0) continue;
        // std::cout << i << std::endl;
        Array::array_xxc c_data = Array::array_xxc::Zero(end_ch - start_ch + 2 * npad_wire, m_end_tick - m_start_tick);

        // fill normal order
        for (size_t j = 0; j != m_vec_vec_charge.at(i).size(); j++) {
            c_data(std::get<0>(m_vec_vec_charge.at(i).at(j)) + npad_wire - start_ch,
                   std::get<1>(m_vec_vec_charge.at(i).at(j)) - m_start_tick) +=
                std::get<2>(m_vec_vec_charge.at(i).at(j));
        }
        //    std::cout << i << " " << m_vec_vec_charge.at(i).size() << std::endl;
        m_vec_vec_charge.at(i).clear();
        m_vec_vec_charge.at(i).shrink_to_fit();

        // fill reverse order
        int ii = num_double * 2 - i;
        for (size_t j = 0; j != m_vec_vec_charge.at(ii).size(); j++) {
            c_data(end_ch + npad_wire - 1 - std::get<0>(m_vec_vec_charge.at(ii).at(j)),
                   std::get<1>(m_vec_vec_charge.at(ii).at(j)) - m_start_tick) +=
                std::complex<float>(0, std::get<2>(m_vec_vec_charge.at(ii).at(j)));
        }
        //    std::cout << ii << " " << m_vec_vec_charge.at(ii).size() << std::endl;
        m_vec_vec_charge.at(ii).clear();
        m_vec_vec_charge.at(ii).shrink_to_fit();

        // Do FFT on time
        // Do FFT on wire
        c_data = fwd(m_dft, c_data);

        // std::cout << i << std::endl;
        {
            Array::array_xxc resp_f_w =
                Array::array_xxc::Zero(end_ch - start_ch + 2 * npad_wire, m_end_tick - m_start_tick);
            {
                Waveform::compseq_t rs1 = m_vec_map_resp.at(i)[0]->spectrum();
                // do a inverse FFT
                Waveform::realseq_t rs1_t = inv_c2r(m_dft, rs1);
                // pick the first xxx ticks
                Waveform::realseq_t rs1_reduced(m_end_tick - m_start_tick, 0);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    if (icol >= int(rs1_t.size())) break;
                    rs1_reduced.at(icol) = rs1_t[icol];
                }
                // do a FFT
                rs1 = fwd_r2c(m_dft, rs1_reduced);

                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    resp_f_w(0, icol) = rs1[icol];
                }
            }

            for (int irow = 0; irow != m_num_pad_wire; irow++) {
                Waveform::compseq_t rs1 = m_vec_map_resp.at(i)[irow + 1]->spectrum();
                Waveform::realseq_t rs1_t = inv_c2r(m_dft, rs1);
                Waveform::realseq_t rs1_reduced(m_end_tick - m_start_tick, 0);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    if (icol >= int(rs1_t.size())) break;
                    rs1_reduced.at(icol) = rs1_t[icol];
                }
                rs1 = fwd_r2c(m_dft, rs1_reduced);
                Waveform::compseq_t rs2 = m_vec_map_resp.at(i)[-irow - 1]->spectrum();
                Waveform::realseq_t rs2_t = inv_c2r(m_dft, rs2);
                Waveform::realseq_t rs2_reduced(m_end_tick - m_start_tick, 0);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    if (icol >= int(rs2_t.size())) break;
                    rs2_reduced.at(icol) = rs2_t[icol];
                }
                rs2 = fwd_r2c(m_dft, rs2_reduced);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    resp_f_w(irow + 1, icol) = rs1[icol];
                    resp_f_w(end_ch - start_ch - 1 - irow + 2 * npad_wire, icol) = rs2[icol];
                }
            }
            // std::cout << i << std::endl;

            // Do FFT on wire for response // slight larger
            // Now becomes the f and f in both time and wire domain ...
            resp_f_w = fwd(m_dft, resp_f_w, 0);

            // multiply them together
            c_data = c_data * resp_f_w;
        }

        // Do inverse FFT on wire
        c_data = inv(m_dft, c_data, 0);

        // Add to wire result in frequency
        acc_data_f_w += c_data;
    }

    // std::cout << "ABC : " << std::endl;

    // central region ...
    {
        int i = num_double;
        // fill response array in frequency domain

        Array::array_xxc data_f_w;
        {
            Array::array_xxf data_t_w =
                Array::array_xxf::Zero(end_ch - start_ch + 2 * npad_wire, m_end_tick - m_start_tick);
            // fill charge array in time-wire domain // slightly larger
            for (size_t j = 0; j != m_vec_vec_charge.at(i).size(); j++) {
                data_t_w(std::get<0>(m_vec_vec_charge.at(i).at(j)) + npad_wire - start_ch,
                         std::get<1>(m_vec_vec_charge.at(i).at(j)) - m_start_tick) +=
                    std::get<2>(m_vec_vec_charge.at(i).at(j));
                // std::cout << std::get<1>(m_vec_vec_charge.at(i).at(j)) << std::endl;
            }
            //      std::cout << i << " " << m_vec_vec_charge.at(i).size() << std::endl;
            m_vec_vec_charge.at(i).clear();
            m_vec_vec_charge.at(i).shrink_to_fit();

            // Do FFT on time
            data_f_w = fwd_r2c(m_dft, data_t_w, 1);
            // Do FFT on wire
            data_f_w = fwd(m_dft, data_f_w, 0);
        }

        {
            Array::array_xxc resp_f_w =
                Array::array_xxc::Zero(end_ch - start_ch + 2 * npad_wire, m_end_tick - m_start_tick);

            {
                Waveform::compseq_t rs1 = m_vec_map_resp.at(i)[0]->spectrum();

                // do a inverse FFT
                Waveform::realseq_t rs1_t = inv_c2r(m_dft, rs1);
                // pick the first xxx ticks
                Waveform::realseq_t rs1_reduced(m_end_tick - m_start_tick, 0);
                // std::cout << rs1.size() << " " << nsamples << " " << m_end_tick << " " <<  m_start_tick << std::endl;
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    if (icol >= int(rs1_t.size())) break;
                    rs1_reduced.at(icol) = rs1_t[icol];
                    //  std::cout << icol << " " << rs1_t[icol] << std::endl;
                }
                // do a FFT
                rs1 = fwd_r2c(m_dft, rs1_reduced);

                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    //   std::cout << icol << " " << rs1[icol] << " " << temp_resp_f_w(0,icol) << std::endl;
                    resp_f_w(0, icol) = rs1[icol];
                }
            }
            for (int irow = 0; irow != m_num_pad_wire; irow++) {
                Waveform::compseq_t rs1 = m_vec_map_resp.at(i)[irow + 1]->spectrum();
                Waveform::realseq_t rs1_t = inv_c2r(m_dft, rs1);
                Waveform::realseq_t rs1_reduced(m_end_tick - m_start_tick, 0);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    if (icol >= int(rs1_t.size())) break;
                    rs1_reduced.at(icol) = rs1_t[icol];
                }
                rs1 = fwd_r2c(m_dft, rs1_reduced);
                Waveform::compseq_t rs2 = m_vec_map_resp.at(i)[-irow - 1]->spectrum();
                Waveform::realseq_t rs2_t = inv_c2r(m_dft, rs2);
                Waveform::realseq_t rs2_reduced(m_end_tick - m_start_tick, 0);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    if (icol >= int(rs2_t.size())) break;
                    rs2_reduced.at(icol) = rs2_t[icol];
                }
                rs2 = fwd_r2c(m_dft, rs2_reduced);
                for (int icol = 0; icol != m_end_tick - m_start_tick; icol++) {
                    resp_f_w(irow + 1, icol) = rs1[icol];
                    resp_f_w(end_ch - start_ch - 1 - irow + 2 * npad_wire, icol) = rs2[icol];
                }
            }
            // Do FFT on wire for response // slight larger
            // Now becomes the f and f in both time and wire domain ...
            resp_f_w = fwd(m_dft, resp_f_w, 0);
            // multiply them together
            data_f_w = data_f_w * resp_f_w;
        }

        // Do inverse FFT on wire
        data_f_w = inv(m_dft, data_f_w, 0);

        // Add to wire result in frequency
        acc_data_f_w += data_f_w;
    }

    acc_data_f_w = inv(m_dft, acc_data_f_w, 1); 
    Array::array_xxf real_m_decon_data = acc_data_f_w.real();
    Array::array_xxf img_m_decon_data = acc_data_f_w.imag().colwise().reverse();
    m_decon_data = real_m_decon_data + img_m_decon_data;

}  // constructor

Gen::ImpactTransform::~ImpactTransform() {}

Waveform::realseq_t Gen::ImpactTransform::waveform(int iwire) const
{
    const int nsamples = m_bd.tbins().nbins();
    if (iwire < m_start_ch || iwire >= m_end_ch) {
        return Waveform::realseq_t(nsamples, 0.0);
    }
    else {
        Waveform::realseq_t wf(nsamples, 0.0);
        for (int i = 0; i != nsamples; i++) {
            if (i >= m_start_tick && i < m_end_tick) {
                wf.at(i) = m_decon_data(iwire - m_start_ch, i - m_start_tick);
            }
            // std::cout << m_decon_data(iwire-m_start_ch,i-m_start_tick) << std::endl;
        }

        if (m_pir->closest(0)->long_aux_waveform().size() > 0) {
            // now convolute with the long-range response ...
            const size_t nlength = fft_best_length(nsamples + m_pir->closest(0)->long_aux_waveform_pad());

            //   std::cout << nlength << " " << nsamples + m_pir->closest(0)->long_aux_waveform_pad() << std::endl;

            wf.resize(nlength, 0);
            Waveform::realseq_t long_resp = m_pir->closest(0)->long_aux_waveform();
            long_resp.resize(nlength, 0);
            Waveform::compseq_t spec = fwd_r2c(m_dft, wf);
            Waveform::compseq_t long_spec = fwd_r2c(m_dft, long_resp);
            for (size_t i = 0; i != nlength; i++) {
                spec.at(i) *= long_spec.at(i);
            }
            wf = inv_c2r(m_dft, spec);
            wf.resize(nsamples, 0);
        }

        return wf;
    }
}
