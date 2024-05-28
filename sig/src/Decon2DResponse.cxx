#include "WireCellSig/Decon2DResponse.h"

#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/Util.h"
#include "WireCellAux/TensUtil.h"
#include "WireCellAux/DftTools.h"

#include "WireCellIface/ITensorSet.h"
#include "WireCellIface/IFilterWaveform.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Response.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellUtil/Exceptions.h"


WIRECELL_FACTORY(Decon2DResponse, WireCell::Sig::Decon2DResponse, WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::DftTools::fwd;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv;
using WireCell::Aux::DftTools::inv_c2r;


Sig::Decon2DResponse::Decon2DResponse()
  : log(Log::logger("sig"))
{
}

Configuration Sig::Decon2DResponse::default_configuration() const
{
    Configuration cfg;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    return cfg;
}

void Sig::Decon2DResponse::configure(const WireCell::Configuration &cfg)
{
    m_cfg = cfg;

    m_anode = Factory::find_tn<IAnodePlane>(get<std::string>(m_cfg, "anode", "Anode"));
    if (!m_anode) {
        THROW(ValueError() << errmsg{"Sig::Decon2DResponse::configure !m_anode"});
    }

    auto cresp_n = get<std::string>(m_cfg, "per_chan_resp", "PerChannelResponse");
    if (!cresp_n.empty()) {
        m_cresp = Factory::find_tn<IChannelResponse>(cresp_n);
        if (!m_cresp) {
            THROW(ValueError() << errmsg{"Sig::Decon2DResponse::configure !m_cresp"});
        }
    }
    else {
        m_cresp = nullptr;
        log->warn("Sig::Decon2DResponse::configure No PerChannelResponse");
    }

    m_fresp = Factory::find_tn<IFieldResponse>(get<std::string>(m_cfg, "field_response", "FieldResponse"));
    if (!m_fresp) {
        THROW(ValueError() << errmsg{"Sig::Decon2DResponse::configure !m_fresp"});
    }

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
}

// unused
#if 0
namespace {
    std::string dump(const ITensorSet::pointer &itens)
    {
        std::stringstream ss;
        ss << "ITensorSet: ";
        ss << itens->ident() << ", " << itens->metadata();
        for (auto iten : *itens->tensors()) {
            ss << "shape: [";
            for (auto l : iten->shape()) {
                ss << l << " ";
            }
            ss << "]\n";
        }
        return ss.str();
    }
}  // namespace
#endif

std::vector<Waveform::realseq_t> Sig::Decon2DResponse::init_overall_response(const ITensorSet::pointer &in) const
{
    // TODO: protection?
    auto set_md = in->metadata();
    auto tag = get<std::string>(m_cfg, "tag", "trace_tag");
    auto iwf_ten = Aux::get_tens(in, tag, "waveform");
    auto wf_shape = iwf_ten->shape();
    // size_t nchans = wf_shape[0];
    size_t nticks = wf_shape[1];

    // bins
    double m_period = get(set_md, "tick", 0.5 * units::microsecond);
    // int m_fft_nwires = nchans;
    // int m_pad_nwires = 0;
    int m_fft_nticks = nticks;
    // int m_pad_nticks = 0;

    // average overall responses
    std::vector<Waveform::realseq_t> m_overall_resp;

    // gain, shaping time, other applification factors
    double m_gain = get<double>(m_cfg, "gain", 14.0 * units::mV / units::fC);
    double m_shaping_time = get<double>(m_cfg, "shaping_time", 2.2 * units::microsecond);
    double m_inter_gain = get<double>(m_cfg, "inter_gain", 1.2);
    double m_ADC_mV = get<double>(m_cfg, "ADC_mV", 4096 / (2000. * units::mV));

    // time offset
    double m_fine_time_offset = get<double>(m_cfg, "ftoffset", 0.0 * units::microsecond);
    // double m_coarse_time_offset = get<double>(m_cfg, "ctoffset", -8.0 * units::microsecond);

    // Get full, "fine-grained" field responses defined at impact positions.
    Response::Schema::FieldResponse fr = m_fresp->field_response();

    // Make a new data set which is the average FR
    Response::Schema::FieldResponse fravg = Response::wire_region_average(fr);

    // since we only do FFT along time, no need to change dimension for wire ...
    const size_t fine_nticks = fft_best_length(fravg.planes[0].paths[0].current.size());
    int fine_nwires = fravg.planes[0].paths.size();

    WireCell::Waveform::compseq_t elec;
    WireCell::Binning tbins(fine_nticks, 0, fine_nticks * fravg.period);
    Response::ColdElec ce(m_gain, m_shaping_time);
    auto ewave = ce.generate(tbins);
    Waveform::scale(ewave, m_inter_gain * m_ADC_mV * (-1));
    elec = fwd_r2c(m_dft, ewave);

    std::complex<float> fine_period(fravg.period, 0);

    Waveform::realseq_t wfs(m_fft_nticks);
    Waveform::realseq_t ctbins(m_fft_nticks);
    for (int i = 0; i != m_fft_nticks; i++) {
        ctbins.at(i) = i * m_period;
    }

    Waveform::realseq_t ftbins(fine_nticks);
    for (size_t i = 0; i != fine_nticks; i++) {
        ftbins.at(i) = i * fravg.period;
    }

    // clear the overall response
    m_overall_resp.clear();

    // Convert each average FR to a 2D array
    int iplane = get<int>(m_cfg, "iplane", 0);
    auto arr = Response::as_array(fravg.planes[iplane], fine_nwires, fine_nticks);

    // do FFT for response ...
    Array::array_xxc c_data = fwd_r2c(m_dft, arr, 1);

    int nrows = c_data.rows();
    int ncols = c_data.cols();

    for (int irow = 0; irow < nrows; ++irow) {
        for (int icol = 0; icol < ncols; ++icol) {
            c_data(irow, icol) = c_data(irow, icol) * elec.at(icol) * fine_period;
        }
    }

    arr = inv_c2r(m_dft, c_data, 1);

    // figure out how to do fine ... shift (good ...)
    int fine_time_shift = m_fine_time_offset / fravg.period;
    if (fine_time_shift > 0) {
        Array::array_xxf arr1(nrows, ncols - fine_time_shift);
        arr1 = arr.block(0, 0, nrows, ncols - fine_time_shift);
        Array::array_xxf arr2(nrows, fine_time_shift);
        arr2 = arr.block(0, ncols - fine_time_shift, nrows, fine_time_shift);
        arr.block(0, 0, nrows, fine_time_shift) = arr2;
        arr.block(0, fine_time_shift, nrows, ncols - fine_time_shift) = arr1;
    }

    // redigitize ...
    for (int irow = 0; irow < fine_nwires; ++irow) {
        size_t fcount = 1;
        for (int i = 0; i != m_fft_nticks; i++) {
            double ctime = ctbins.at(i);
            if (fcount < fine_nticks)
                while (ctime > ftbins.at(fcount)) {
                    fcount++;
                    if (fcount >= fine_nticks) break;
                }

            if (fcount < fine_nticks) {
                wfs.at(i) = ((ctime - ftbins.at(fcount - 1)) / fravg.period * arr(irow, fcount - 1) +
                             (ftbins.at(fcount) - ctime) / fravg.period * arr(irow, fcount));  // / (-1);
            }
            else {
                wfs.at(i) = 0;
            }
        }
        m_overall_resp.push_back(wfs);
    }  // loop inside wire ...

    return m_overall_resp;
}

bool Sig::Decon2DResponse::operator()(const ITensorSet::pointer &in, ITensorSet::pointer &out)
{
    log->debug("Decon2DResponse: start");

    if (!in) {
        out = nullptr;
        return true;
    }
    auto set_md = in->metadata();

    auto tag = get<std::string>(m_cfg, "tag", "trace_tag");

    auto iwf_ten = Aux::get_tens(in, tag, "waveform");
    if (!iwf_ten) {
        log->error("Tensor->Frame: failed to get waveform tensor for tag {}", tag);
        THROW(ValueError() << errmsg{"Sig::Decon2DResponse::operator() !iwf_ten"});
    }
    auto wf_shape = iwf_ten->shape();
    size_t nchans = wf_shape[0];
    size_t nticks = wf_shape[1];
    log->debug("iwf_ten->shape: {} {}", nchans, nticks);

    auto ch_ten = Aux::get_tens(in, tag, "channels");
    if (!ch_ten) {
        log->error("Tensor->Frame: failed to get channels tensor for tag {}", tag);
        THROW(ValueError() << errmsg{"Sig::Decon2DResponse::operator() !ch_ten"});
    }
    Eigen::Map<const Eigen::ArrayXi> ch_arr((const int *) ch_ten->data(), nchans);

    bool have_summary = true;
    auto sum_ten = Aux::get_tens(in, tag, "summary");
    if (!sum_ten) {
        log->debug("Tensor->Frame: no optional summary tensor for tag {}", tag);
        have_summary = false;
    }

    int iplane = get<int>(m_cfg, "iplane", 0);

    const std::vector<std::string> filter_names{"Wire_ind", "Wire_ind", "Wire_col"};

    // bins
    double m_period = get(set_md, "tick", 0.5 * units::microsecond);
    // int m_fft_nwires = nchans;
    // int m_pad_nwires = 0;
    int m_fft_nticks = nticks;
    // int m_pad_nticks = 0;

    // time offset
    // double m_fine_time_offset = get<double>(m_cfg, "ftoffset", 0.0 * units::microsecond);
    double m_coarse_time_offset = get<double>(m_cfg, "ctoffset", -8.0 * units::microsecond);
    auto fr = m_fresp->field_response();
    double m_intrinsic_time_offset = fr.origin / fr.speed;

    // gain, shaping time, other applification factors
    double m_gain = get<double>(m_cfg, "gain", 14.0 * units::mV / units::fC);
    double m_shaping_time = get<double>(m_cfg, "shaping_time", 2.2 * units::microsecond);
    // double m_inter_gain = get<double>(m_cfg, "inter_gain", 1.2);
    // double m_ADC_mV = get<double>(m_cfg, "ADC_mV", 4096 / (2000. * units::mV));

    auto m_overall_resp = init_overall_response(in);
    log->debug("m_overall_resp: {}", m_overall_resp.size());

    // calculated the wire shift ...
    int m_wire_shift = (int(m_overall_resp.size()) - 1) / 2;

    // TensorSet to Eigen
    WireCell::Array::array_xxf r_data = Aux::itensor_to_eigen_array<float>(iwf_ten);
    log->debug("r_data: {} {}", r_data.rows(), r_data.cols());

    // first round of FFT on time
    WireCell::Array::array_xxc c_data = fwd_r2c(m_dft, r_data, 1);

    if (m_cresp) {
        log->debug("Decon2DResponse: applying ch-by-ch electronics response correction");
        auto cr_bins = m_cresp->channel_response_binning();
        if (cr_bins.binsize() != m_period) {
            log->critical("Decon2DResponse::decon_2D_init: channel response size mismatch");
            THROW(ValueError() << errmsg{"Decon2DResponse::decon_2D_init: channel response size mismatch"});
        }
        WireCell::Binning tbins(m_fft_nticks, cr_bins.min(), cr_bins.min() + m_fft_nticks * m_period);
        Response::ColdElec ce(m_gain, m_shaping_time);

        const auto ewave = ce.generate(tbins);
        const WireCell::Waveform::compseq_t elec = fwd_r2c(m_dft, ewave);

        for (int irow = 0; irow != c_data.rows(); irow++) {
            Waveform::realseq_t tch_resp = m_cresp->channel_response(ch_arr[irow]);
            tch_resp.resize(m_fft_nticks, 0);
            const WireCell::Waveform::compseq_t ch_elec = fwd_r2c(m_dft, tch_resp);

            // FIXME figure this out
            // const int irow = och.wire + m_pad_nwires;
            for (int icol = 0; icol != c_data.cols(); icol++) {
                const auto four = ch_elec.at(icol);
                if (std::abs(four) != 0) {
                    c_data(irow, icol) *= elec.at(icol) / four;
                }
                else {
                    c_data(irow, icol) = 0;
                }
            }
        }
    }
    log->trace("TRACE {}", __LINE__);

    // second round of FFT on wire
    c_data = fwd(m_dft, c_data, 0);

    // response part ...
    Array::array_xxf r_resp = Array::array_xxf::Zero(r_data.rows(), m_fft_nticks);
    for (size_t i = 0; i != m_overall_resp.size(); i++) {
        for (int j = 0; j != m_fft_nticks; j++) {
            r_resp(i, j) = m_overall_resp.at(i).at(j);
        }
    }
    log->trace("TRACE {}", __LINE__);

    // do first round FFT on the resposne on time
    Array::array_xxc c_resp = fwd_r2c(m_dft, r_resp, 1);
    // do second round FFT on the response on wire
    c_resp = fwd(m_dft, c_resp, 0);

    // make ratio to the response and apply wire filter
    c_data = c_data / c_resp;
    log->trace("TRACE {}", __LINE__);

    // apply software filter on wire
    Waveform::realseq_t wire_filter_wf;
    auto ncr1 = Factory::find<IFilterWaveform>("HfFilter", filter_names[iplane]);
    wire_filter_wf = ncr1->filter_waveform(c_data.rows());
    for (int irow = 0; irow < c_data.rows(); ++irow) {
        for (int icol = 0; icol < c_data.cols(); ++icol) {
            float val = abs(c_data(irow, icol));
            if (std::isnan(val)) {
                c_data(irow, icol) = -0.0;
            }
            if (std::isinf(val)) {
                c_data(irow, icol) = 0.0;
            }
            c_data(irow, icol) *= wire_filter_wf.at(irow);
        }
    }
    log->trace("TRACE {}", __LINE__);

    // do the first round of inverse FFT on wire
    c_data = inv(m_dft, c_data, 0);

    // do the second round of inverse FFT on time
    r_data = inv_c2r(m_dft, c_data, 1);

    // do the shift in wire
    const int nrows = r_data.rows();
    const int ncols = r_data.cols();
    {
        Array::array_xxf arr1(m_wire_shift, ncols);
        arr1 = r_data.block(nrows - m_wire_shift, 0, m_wire_shift, ncols);
        Array::array_xxf arr2(nrows - m_wire_shift, ncols);
        arr2 = r_data.block(0, 0, nrows - m_wire_shift, ncols);
        r_data.block(0, 0, m_wire_shift, ncols) = arr1;
        r_data.block(m_wire_shift, 0, nrows - m_wire_shift, ncols) = arr2;
    }

    // do the shift in time
    int time_shift = (m_coarse_time_offset + m_intrinsic_time_offset) / m_period;
    if (time_shift > 0) {
        Array::array_xxf arr1(nrows, ncols - time_shift);
        arr1 = r_data.block(0, 0, nrows, ncols - time_shift);
        Array::array_xxf arr2(nrows, time_shift);
        arr2 = r_data.block(0, ncols - time_shift, nrows, time_shift);
        r_data.block(0, 0, nrows, time_shift) = arr2;
        r_data.block(0, time_shift, nrows, ncols - time_shift) = arr1;
    }
    c_data = fwd_r2c(m_dft, r_data, 1);

    log->trace("TRACE {}", __LINE__);

    // Eigen to TensorSet
    auto owf_ten = Aux::eigen_array_to_simple_tensor<std::complex<float> >(c_data);

    // Debug: inverte fft and save out
    // auto owf_ten = Aux::eigen_array_to_simple_tensor<float>(r_data);

    // save back to ITensorSet
    ITensor::vector *itv = new ITensor::vector;
    auto iwf_md = iwf_ten->metadata();
    auto &owf_md = owf_ten->metadata();
    owf_md["tag"] = tag;
    owf_md["pad"] = get<float>(iwf_md, "pad", 0.0);
    owf_md["tbin"] = get<float>(iwf_md, "tbin", 0.0);
    owf_md["type"] = "waveform";
    itv->push_back(ITensor::pointer(owf_ten));

    // FIXME need to change ch_ten tag to outtag
    itv->push_back(ch_ten);
    if (have_summary) {
        itv->push_back(sum_ten);
    }

    Configuration oset_md(in->metadata());
    oset_md["tags"] = Json::arrayValue;
    oset_md["tags"].append(tag);

    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), oset_md, ITensor::shared_vector(itv));

    log->debug("Decon2DResponse: end");

    return true;
}
