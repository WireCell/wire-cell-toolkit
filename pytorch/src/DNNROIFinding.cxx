
#include "WireCellPytorch/DNNROIFinding.h"
#include "WireCellPytorch/Util.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/ITrace.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/PlaneTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/TimeKeeper.h"

#include <string>
#include <vector>

// Uncomment to enable local writing of HDF5 file for debugging
// purposes.  A temporary HDF5 dependency must be added in
// wscript_build.  Do not commit with this.

// #define DNNROI_HDF5_DEBUG

// FIXME: consider changing this to instead write out numpy.tar.bz2
// streams which does not require HDF5 dependency and could then
// become a dynamic option.

#ifdef DNNROI_HDF5_DEBUG
#include <h5cpp/all>
#endif

/// macro to register name - concrete pair in the NamedFactory
/// @param NAME - used to configure node in JSON/Jsonnet
/// @parame CONCRETE - C++ concrete type
/// @parame ... - interfaces
WIRECELL_FACTORY(DNNROIFinding, WireCell::Pytorch::DNNROIFinding, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;

Pytorch::DNNROIFinding::DNNROIFinding()
    : Aux::Logger("DNNROIFinding", "torch")
    , m_save_count(0)
{
}

Pytorch::DNNROIFinding::~DNNROIFinding() {}

void Pytorch::DNNROIFinding::configure(const WireCell::Configuration& cfg)
{
    m_cfg.anode = get(cfg, "anode", m_cfg.anode);
    m_cfg.plane = get(cfg, "plane", m_cfg.plane);

    auto anode = Factory::find_tn<IAnodePlane>(m_cfg.anode);
    auto ichans = Aux::plane_channels(anode, m_cfg.plane);
    m_chset.clear();
    m_chlist.clear();
    for (const auto& ichan : ichans) {
        auto chid = ichan->ident();
        m_chset.insert(chid);
        m_chlist.push_back(chid);
    }
    {
        auto chmm = std::minmax_element(m_chset.begin(), m_chset.end());
        log->debug("anode={} plane={} nchans={} in=[{},{}]",
                   m_cfg.anode, m_cfg.plane, m_chset.size(),
                   *chmm.first, *chmm.second);
    }
    m_cfg.sort_chanids = get(cfg, "sort_chanids", m_cfg.sort_chanids);
    if (m_cfg.sort_chanids) {
        std::sort(m_chlist.begin(), m_chlist.end());
    }

    m_cfg.input_scale = get(cfg, "input_scale", m_cfg.input_scale);
    m_cfg.input_offset = get(cfg, "input_offset", m_cfg.input_offset);
    m_cfg.output_scale = get(cfg, "output_scale", m_cfg.output_scale);
    m_cfg.output_offset = get(cfg, "output_offset", m_cfg.output_offset);
    if (m_cfg.output_scale != 1.0) {
        log->debug("using output scale: {}", m_cfg.output_scale);
    }
        
    m_cfg.tick0 = get(cfg, "tick0", m_cfg.tick0);
    m_cfg.nticks = get(cfg, "nticks", m_cfg.nticks);
    m_cfg.mask_thresh = get(cfg, "mask_thresh", m_cfg.mask_thresh);
    m_cfg.forward = get(cfg, "forward", m_cfg.forward);
    m_cfg.intags.clear();
    for (auto one : cfg["intags"]) {
        m_cfg.intags.push_back(one.asString());
    }
    m_cfg.summary_tag = get(cfg, "summary_tag", m_cfg.summary_tag);
    m_cfg.tick_per_slice = get(cfg, "tick_per_slice", m_cfg.tick_per_slice);
    m_cfg.decon_charge_tag = get(cfg, "decon_charge_tag", m_cfg.decon_charge_tag);
    m_cfg.outtag = get(cfg, "outtag", m_cfg.outtag);
    m_cfg.debugfile = get(cfg, "debugfile", m_cfg.debugfile);
    m_cfg.nchunks = get(cfg, "nchunks", m_cfg.nchunks);
    m_cfg.save_negative_charge = get(cfg, "save_negative_charge", m_cfg.save_negative_charge);

    m_nrows = m_chlist.size();
    m_ncols = m_cfg.nticks;

    if (m_cfg.intags.empty()) {
        log->critical("intags is empty");
        THROW(ValueError() << errmsg{"intags is empty"});
    }
    if (m_cfg.decon_charge_tag.empty()) {
        log->critical("decon_charge_tag is empty");
        THROW(ValueError() << errmsg{"decon_charge_tag is empty"});
    }
    if (!cfg["cbeg"].isNull() or !cfg["cend"].isNull()) {
        log->warn("the cbeg/cend option is not supported, use 'plane' number");
    }
    if (m_nrows == 0) {
        log->critical("empty channel ID list, configure 'anode' and 'plane' parameter");
        THROW(ValueError() << errmsg{"empty channel ID list"});
    }

    m_forward = Factory::find_tn<ITensorForward>(m_cfg.forward);

    // We output dense traces always of same size so make the trace
    // indices once for all later reuse.
    m_trace_indices.resize(m_nrows);
    std::iota(m_trace_indices.begin(), m_trace_indices.end(), 0);


#ifdef DNNROI_HDF5_DEBUG
    {
        if (m_cfg.debugfile.empty()) {
            log->critical("Must provide output evalfile to DNNROIFinding");
            THROW(ValueError() << errmsg{"Must provide output evalfile to DNNROIFinding"});
        }
        h5::create(m_cfg.debugfile, H5F_ACC_TRUNC);
    }
#endif

}

WireCell::Configuration Pytorch::DNNROIFinding::default_configuration() const
{
    Configuration cfg;

    // see coments under DNNROIFindingCfg in header file
    cfg["anode"] = m_cfg.anode;
    cfg["plane"] = m_cfg.plane;
    cfg["sort_chanids"] = m_cfg.sort_chanids;

    cfg["input_scale"] = m_cfg.input_scale;
    cfg["input_offset"] = m_cfg.input_offset;
    cfg["output_scale"] = m_cfg.output_scale;
    cfg["output_offset"] = m_cfg.output_offset;
    cfg["tick0"] = m_cfg.tick0;
    cfg["nticks"] = m_cfg.nticks;
    cfg["mask_thresh"] = m_cfg.mask_thresh;
    cfg["forward"] = m_cfg.forward;
    cfg["intags"] = Json::arrayValue;
    for (const auto& one : m_cfg.intags) {
        cfg["intags"].append(one);
    }
    cfg["tick_per_slice"] = m_cfg.tick_per_slice;
    cfg["decon_charge_tag"] = m_cfg.decon_charge_tag;
    cfg["outtag"] = m_cfg.outtag;
    cfg["debugfile"] = m_cfg.debugfile;
    cfg["nchunks"] = m_cfg.nchunks;
    cfg["save_negative_charge"] = m_cfg.save_negative_charge;
    return cfg;
}


ITrace::vector Pytorch::DNNROIFinding::select(ITrace::vector traces)
{
    auto end = std::remove_if(traces.begin(), traces.end(),
                              [&](const ITrace::pointer& t) {
                                  return m_chset.find(t->channel()) == m_chset.end();
                              });
    traces.resize(end - traces.begin());
    return traces;
}


Array::array_xxf Pytorch::DNNROIFinding::traces_to_eigen(ITrace::vector traces)
{
    Array::array_xxf arr = Array::array_xxf::Zero(m_nrows, m_ncols);
    traces = select(traces);
    if (traces.empty()) {
        return arr;
    }
    Aux::fill(arr, traces, m_chlist.begin(), m_chlist.end(), m_cfg.tick0);
    return arr;
}

IFrame::trace_summary_t Pytorch::DNNROIFinding::get_summary_e(const IFrame::pointer& inframe, const std::string &tag) const {
    IFrame::trace_summary_t summary_e(m_nrows, 0.0);
    std::unordered_map<int, size_t> ch2row;
    for (size_t i = 0; i < m_chlist.size(); ++i) {
        ch2row[m_chlist[i]] = i;
    }
    auto traces = Aux::tagged_traces(inframe, tag);
    auto summary = inframe->trace_summary(tag);
    for (size_t i = 0; i < traces.size(); ++i) {
        const auto& tr = traces[i];
        const auto& ch = tr->channel();
        const auto& row = ch2row.find(ch);
        if (row != ch2row.end()) {
            summary_e[row->second] = summary[i];
        }
    }
    return summary_e;
}

ITrace::shared_vector Pytorch::DNNROIFinding::eigen_to_traces(const Array::array_xxf& arr, bool set_negative_to_zero)
{
    ITrace::vector traces;
    ITrace::ChargeSequence charge(m_ncols, 0.0);
    for (size_t irow = 0; irow < m_nrows; ++irow) {
        auto wave = arr.row(irow);
        for (size_t icol=0; icol<m_ncols; ++icol) {
            charge[icol] = wave(icol);
            if (set_negative_to_zero) {
                charge[icol] = charge[icol] < 0 ? 0 : charge[icol];
            }
        }
        const auto ch = m_chlist[irow];
        traces.push_back(std::make_shared<Aux::SimpleTrace>(ch, 0, charge));
    }
    return std::make_shared<ITrace::vector>(traces.begin(), traces.end());
}


bool Pytorch::DNNROIFinding::operator()(const IFrame::pointer& inframe, IFrame::pointer& outframe)
{
    outframe = inframe;
    if (!inframe) {
        log->debug("EOS at call={}", m_save_count);
        outframe = nullptr;
        return true;
    }

    TimeKeeper tk(fmt::format("call={}", m_save_count));

    // frame to eigen
    std::vector<Array::array_xxf> ch_eigen;
    for (auto tag : m_cfg.intags) {
        auto traces = Aux::tagged_traces(inframe, tag);
        auto arr = traces_to_eigen(traces);
        if (arr.sum() == 0.0) {
            log->warn("call={} no traces for input tag {}, using zeros", m_save_count, tag);
        }
        else { 
            log->debug("call={} tag={} ntraces={}", m_save_count, tag, traces.size());
        }
        arr = arr * m_cfg.input_scale + m_cfg.input_offset;
        ch_eigen.push_back(Array::downsample(arr, m_cfg.tick_per_slice, 1));
    }

    // eigen to tensor
    std::vector<torch::Tensor> ch;
    for (unsigned int i = 0; i < ch_eigen.size(); ++i) {
        ch.push_back(torch::from_blob(ch_eigen[i].data(), {ch_eigen[i].cols(), ch_eigen[i].rows()}));
    }
    // ret: {ntags, nticks, nchannels}
    auto img = torch::stack(ch, 0);
    // ret: {1, ntags, nchannels, nticks}
    auto batch = torch::stack({torch::transpose(img, 1, 2)}, 0);

    auto chunks = batch.chunk(m_cfg.nchunks, 2);
    std::vector<torch::Tensor> outputs;

    log->debug(tk(fmt::format("call={} calling model \"{}\" with {} chunks ",
                              m_save_count, m_cfg.forward, m_cfg.nchunks)));
    for (auto chunk : chunks) {
        std::vector<torch::IValue> itens {(m_cfg.nchunks > 1) ? chunk.clone() : chunk};
        auto iitens = Pytorch::to_itensor(itens);
        auto oitens = m_forward->forward(iitens);
        torch::Tensor ochunk = Pytorch::from_itensor({oitens}).front().toTensor().cpu();
        // NOTE: ochunk DOES NOT copy data from oitens!
        // from_blob() in from_itensor() is used to create a tensor from the same memory.
        outputs.push_back(ochunk.clone());
    }
    torch::Tensor output = torch::cat(outputs, 2);
    log->debug(tk(fmt::format("call={} inference done", m_save_count)));

    // tensor to eigen
    Eigen::Map<Eigen::ArrayXXf> out_e(output[0][0].data_ptr<float>(), output.size(3), output.size(2));
    auto mask_e = Array::upsample(out_e, m_cfg.tick_per_slice, 0);

    log->debug(tk(fmt::format("call={} tensor2eigen", m_save_count)));

    // decon charge frame to eigen
    
    Array::array_xxf decon_charge_eigen;
    {
        auto traces = Aux::tagged_traces(inframe, m_cfg.decon_charge_tag);
        decon_charge_eigen = traces_to_eigen(traces);
        if (decon_charge_eigen.sum() == 0.0) {
            log->warn("call={} no traces for input tag {}, using zeros", m_save_count, m_cfg.decon_charge_tag);
        }
        else { 
            log->debug("call={} tag={} ntraces={}", m_save_count, m_cfg.decon_charge_tag, traces.size());
        }
    }

    // apply ROI
    auto sp_charge_T = Array::mask(decon_charge_eigen.transpose(), mask_e, m_cfg.mask_thresh /*0.7*/);
    sp_charge_T = Array::baseline_subtraction(sp_charge_T) * m_cfg.output_scale + m_cfg.output_offset;
    Array::array_xxf sp_charge = sp_charge_T.transpose();

#ifdef DNNROI_HDF5_DEBUG
    // hdf5 eval
    {
        h5::fd_t fd = h5::open(m_cfg.debugevalfile, H5F_ACC_RDWR);
        const unsigned long ncols = mask_e.cols();
        const unsigned long nrows = mask_e.rows();
        // l->info("ncols: {} nrows: {}", ncols, nrows);
        std::string aname = String::format("/%d/frame_%s%d", m_save_count, "dlroi", m_anode->ident());
        h5::write<float>(fd, aname, mask_e.data(), h5::count{ncols, nrows}, h5::chunk{ncols, nrows} | h5::gzip{2});

        aname = String::format("/%d/frame_%s%d", m_save_count, "dlcharge", m_anode->ident());
        h5::write<float>(fd, aname, sp_charge.data(), h5::count{ncols, nrows}, h5::chunk{ncols, nrows} | h5::gzip{2});
    }
    log->debug(tk(fmt::format("call={} h5", m_save_count)));
#endif

    // eigen to frame
    auto traces = eigen_to_traces(sp_charge, !m_cfg.save_negative_charge);
    Aux::SimpleFrame* sframe = new Aux::SimpleFrame(
        inframe->ident(), inframe->time(),
        traces,
        inframe->tick(), inframe->masks());
    sframe->tag_frame("DNNROIFinding");
    log->debug("call={} getting summary for tag {}", m_save_count, m_cfg.summary_tag);
    auto summary = get_summary_e(inframe, m_cfg.summary_tag);
    std::string ss;
    for (const auto& rms : summary) {
        ss += fmt::format("{:.2f} ", rms);
    }
    log->trace("call={} summary: {}", m_save_count, ss);
    sframe->tag_traces(m_cfg.outtag, m_trace_indices, summary);
    outframe = IFrame::pointer(sframe);
    log->debug("call={} output frame: {}", m_save_count, Aux::taginfo(outframe));

    log->debug(tk(fmt::format("call={} finish", m_save_count)));
    ++m_save_count;

    return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
