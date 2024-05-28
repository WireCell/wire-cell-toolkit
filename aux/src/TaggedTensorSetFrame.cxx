#include "WireCellAux/TaggedTensorSetFrame.h"
#include "WireCellAux/TensUtil.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Eigen.h"

WIRECELL_FACTORY(TaggedTensorSetFrame, WireCell::Aux::TaggedTensorSetFrame, WireCell::ITensorSetFrame,
                 WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

WireCell::Configuration Aux::TaggedTensorSetFrame::default_configuration() const
{
    // Each tag will look up a tensor.
    Configuration cfg;
    cfg["tensors"] = Json::arrayValue;
    return cfg;
}

void Aux::TaggedTensorSetFrame::configure(const WireCell::Configuration& config)
{
    auto log = Log::logger("tens");

    m_tags.clear();
    for (auto jten : config["tensors"]) {
        // no tag means empty tag
        auto tag = get<std::string>(jten, "tag", "");
        log->debug("TaggedTensorSetFrame: look for tag: \"{}\"", tag);
        m_tags.insert(tag);
    }
}

bool Aux::TaggedTensorSetFrame::operator()(const input_pointer& in, output_pointer& out)
{
    auto log = Log::logger("tens");
    out = nullptr;
    if (!in) {  // pass on EOS
        log->debug("TaggedTensorSetFrame: EOS");
        return true;
    }

    auto traces = new std::vector<ITrace::pointer>;

    auto set_md = in->metadata();
    auto jtags = set_md["tags"];
    int ident = get(set_md, "ident", 0);

    // Have to make a frame before setting its summaries
    struct SummaryCache {
        std::string tag;
        IFrame::trace_list_t indices;
        IFrame::trace_summary_t summary;
    };
    std::vector<SummaryCache> summaries;
    size_t index_start = 0;

    for (auto jtag : jtags) {
        std::string tag = jtag.asString();
        log->debug("TaggedTensorSetFrame: check input tag: \"{}\"", tag);
        if (m_tags.find(tag) == m_tags.end()) {
            log->debug("Tensor->Frame: skipping unwanted tag {}", tag);
            continue;
        }

        auto wf_ten = Aux::get_tens(in, tag, "waveform");
        if (!wf_ten) {
            log->warn("Tensor->Frame: failed to get waveform tensor for tag {}", tag);
            continue;
        }

        auto ch_ten = Aux::get_tens(in, tag, "channels");
        if (!ch_ten) {
            log->warn("Tensor->Frame: failed to get channels tensor for tag {}", tag);
            continue;
        }

        auto sum_ten = Aux::get_tens(in, tag, "summary");
        if (!sum_ten) {
            log->debug("Tensor->Frame: no optional summary tensor for tag {}", tag);
        }

        auto wf_shape = wf_ten->shape();
        size_t nchans = wf_shape[0];
        size_t nticks = wf_shape[1];
        int tbin = get<int>(wf_ten->metadata(), "tbin", 0);
        log->debug("Tensor->Frame: #{} '{}': {} chans x {} ticks", ident, tag, nchans, nticks);

        Eigen::Map<const Eigen::ArrayXXf> wf_arr((const float*) wf_ten->data(), nchans, nticks);
        Eigen::Map<const Eigen::ArrayXi> ch_arr((const int*) ch_ten->data(), nchans);

        for (size_t ind = 0; ind < nchans; ++ind) {
            SimpleTrace* st = new SimpleTrace(ch_arr[ind], tbin, nticks);
            auto& q = st->charge();
            for (size_t itick = 0; itick < nticks; ++itick) {
                q[itick] = wf_arr(ind, itick);
            }
            // log->debug("Tensor->Frame: #{}: ch: {}, q: {} filled {}", ind, ch_arr[ind], q.size(), nticks);
            traces->push_back(ITrace::pointer(st));
        }

        IFrame::trace_list_t indices;
        for (size_t ind = 0; ind < nchans; ++ind) {
            indices.push_back(index_start + ind);
        }
        index_start += nchans;

        IFrame::trace_summary_t summary;
        if (sum_ten) {
            Eigen::Map<const Eigen::ArrayXd> sum_arr((const double*) sum_ten->data(), nchans);
            summary.resize(nchans);
            for (size_t ind = 0; ind < nchans; ++ind) {
                summary[ind] = sum_arr(ind);
            }
        }

        summaries.emplace_back(SummaryCache{tag, indices, summary});
    }

    SimpleFrame* frame = new SimpleFrame(ident, get(set_md, "time", 0.0), ITrace::shared_vector(traces),
                                         get(set_md, "tick", 0.5 * units::microsecond));

    for (auto& s : summaries) {
        frame->tag_traces(s.tag, s.indices, s.summary);
    }

    out = IFrame::pointer(frame);
    return true;
}

Aux::TaggedTensorSetFrame::TaggedTensorSetFrame() {}
Aux::TaggedTensorSetFrame::~TaggedTensorSetFrame() {}
