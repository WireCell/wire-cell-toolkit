/*! test TaggedFrameTensorSet
 *
 * Also provides an example of how to "hook" an Eigen array to 2D
 * tensor data provided by an ITensor.
 */

#include "WireCellAux/TaggedFrameTensorSet.h"
#include "WireCellAux/TaggedTensorSetFrame.h"
#include "WireCellAux/TensUtil.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Logging.h"

#include "WireCellUtil/Eigen.h"

using namespace WireCell;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

int main()
{
    Log::add_stdout(true, "debug");
    Log::set_level("debug");
    auto log = Log::logger("test");

    const double tick = 0.5 * units::us;
    // "Signal"
    const int signal_frame_ident = 100;
    const double signal_start_time = 6 * units::ms;
    std::vector<float> signal{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 0.0};

    ITrace::vector signal_traces{
        std::make_shared<SimpleTrace>(1, 10, signal), std::make_shared<SimpleTrace>(1, 20, signal),
        std::make_shared<SimpleTrace>(2, 30, signal), std::make_shared<SimpleTrace>(2, 33, signal)};
    auto iframe1 = std::make_shared<SimpleFrame>(signal_frame_ident, signal_start_time, signal_traces, tick);

    ITensorSet::pointer itensorset = nullptr;
    Aux::TaggedFrameTensorSet tfts;
    auto cfg = tfts.default_configuration();
    cfg["tensors"][0]["tag"] = "";
    tfts.configure(cfg);
    bool ok = tfts(iframe1, itensorset);
    Assert(ok);
    Assert(itensorset);

    Assert(itensorset->ident() == signal_frame_ident);
    auto tens = itensorset->tensors();
    Assert(tens);
    Assert(tens->size() == 2);  // waveform+channels

    for (int ind = 0; ind < 2; ++ind) {
        auto ten = tens->at(ind);
        log->debug("tens{} metadata: {}", ind, ten->metadata());
    }

    auto md = itensorset->metadata();
    log->debug("tens set metadata: {}", md);
    Assert(md["time"].asDouble() == signal_start_time);
    Assert(md["tick"].asDouble() == tick);
    auto jtags = md["tags"];
    Assert(jtags.size() == 1);
    Assert(jtags[0].asString() == "");

    for (auto jtag : jtags) {
        std::string tag = jtag.asString();

        auto wft = Aux::get_tens(itensorset, tag, "waveform");
        Assert(wft);

        auto cht = Aux::get_tens(itensorset, tag, "channels");
        Assert(cht);

        auto sumt = Aux::get_tens(itensorset, tag, "summary");
        Assert(sumt == nullptr);

        auto wf_shape = wft->shape();
        Assert(wf_shape.size() == 2);

        const int nchans = wf_shape[0];
        const int nticks = wf_shape[1];
        Assert(nchans == 2);
        log->debug("shape: {} x {}", nchans, nticks);

        auto wf_md = wft->metadata();
        int tbin = wf_md["tbin"].asInt();
        Assert(tbin == 10);

        auto ch_shape = cht->shape();
        Assert(ch_shape.size() == 1);
        Assert(nchans == (int) ch_shape[0]);
        Eigen::Map<Eigen::ArrayXi> ch_arr((int*) cht->data(), nchans);
        Assert(ch_arr[0] == 1);
        Assert(ch_arr[1] == 2);

        Eigen::Map<Eigen::ArrayXXf> wf_arr((float*) wft->data(), nchans, nticks);
        std::stringstream ss;
        ss << "\nchannels:\ntick\t";
        for (int ichan = 0; ichan < nchans; ++ichan) {
            ss << ch_arr[ichan] << "\t";
        }
        ss << "\n";

        for (int itick = 0; itick < nticks; ++itick) {
            ss << tbin + itick;
            for (int ichan = 0; ichan < nchans; ++ichan) {
                ss << "\t" << wf_arr(ichan, itick);
            }
            ss << "\n";
        }
        log->debug(ss.str());
    }

    //
    // Now convert back to an IFrame
    //

    Aux::TaggedTensorSetFrame ttsf;
    cfg = ttsf.default_configuration();
    cfg["tensors"][0]["tag"] = "";
    ttsf.configure(cfg);

    IFrame::pointer iframe2;
    ok = ttsf(itensorset, iframe2);
    Assert(ok);
    Assert(iframe2);
    {
        auto tt = iframe2->trace_tags();
        Assert(tt.size() == 1);
        Assert(tt[0] == "");
    }
    size_t nticks = 0;
    const size_t nchans = 2;
    {
        auto traces = iframe2->traces();
        Assert(traces->size() == nchans);
        Assert(traces->at(0)->channel() == 1);
        nticks = traces->at(0)->charge().size();
        Assert(traces->at(1)->channel() == 2);
        Assert(nticks == traces->at(1)->charge().size());
    }
    {
        auto tt = iframe2->tagged_traces("");
        Assert(tt.size() == nchans);
    }

    {
        auto traces = iframe2->traces();

        std::stringstream ss;
        ss << "\nchannels:\ntick\t";
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
            ss << traces->at(ichan)->channel() << "\t";
        }
        ss << "\n";
        Assert(nticks == 34);
        for (size_t itick = 0; itick < nticks; ++itick) {
            for (size_t ich = 0; ich < nchans; ++ich) {
                const auto& tr = traces->at(ich);
                auto tbin = tr->tbin();
                ss << tbin + itick << ":" << tr->charge()[itick] << " ";
            }
            ss << "\n";
        }
        log->debug(ss.str());
    }
    return 0;
}
