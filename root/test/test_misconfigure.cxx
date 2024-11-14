#include "WireCellAux/DftTools.h"
#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IDFT.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Binning.h"
#include "WireCellUtil/Response.h"
#include "WireCellUtil/NamedFactory.h"

#include <iostream>

#include "MultiPdf.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2F.h"

using namespace std;
using namespace WireCell;
using namespace WireCell::Test;
using namespace WireCell::Waveform;
using namespace WireCell::Aux::DftTools;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

// Make up a Gaussian charge distribution
std::vector<float> blip(Binning bins, double time = 1 * units::ms, double sigma = 3 * units::us, double mag = 10000.0)
{
    int nticks = bins.nbins();
    std::vector<float> charge(nticks, 0);
    for (int ind = 0; ind < nticks; ++ind) {
        const double t = bins.center(ind);
        const double rel = (t - time) / sigma;
        const double val = mag * exp(-0.5 * rel * rel);
        charge[ind] = val;
    }
    return charge;
}

float maximum(const std::vector<float>& wave)
{
    float m = 0;
    for (auto q : wave) {
        if (q > m) m = q;
    }
    return m;
}

TH1F* plot_wave(MultiPdf& pdf, const std::vector<float>& wave, std::string name, std::string title)
{
    TH1F* hist = new TH1F(name.c_str(), title.c_str(), wave.size(), 0, wave.size() + 1);
    int tick = 0;
    for (auto q : wave) {
        hist->Fill(tick + 0.5, q);
        ++tick;
    }
    hist->SetXTitle("tick");
    hist->SetYTitle("amplitude");
    hist->Draw("hist");
    pdf();
    return hist;
}

TH2F* plot_ratio(MultiPdf& pdf, TH2F* f1, TH2F* f2)
{
    TH2F* out = (TH2F*) f1->Clone("ratio");
    out->SetTitle("Ratio - no time shift correction");
    out->Divide(f2);
    out->SetXTitle("ticks");
    out->SetYTitle("channels");
    out->Draw("colz");
    pdf();
    return out;
}

TH2F* plot_frame(MultiPdf& pdf, IFrame::pointer frame, std::string name, double mtick = -1, double mchan = -1);
TH2F* plot_frame(MultiPdf& pdf, IFrame::pointer frame, std::string name, double mtick, double mchan)
{
    auto traces = frame->traces();

    int maxchan = 0, maxtick = 0;
    for (auto trace : (*traces)) {
        auto& charge = trace->charge();
        int tick = trace->tbin() + charge.size();
        if (tick > maxtick) {
            maxtick = tick;
        }
        int chan = trace->channel();
        if (chan > maxchan) {
            maxchan = chan;
        }
        cerr << name << " " << chan << " " << tick << " " << maximum(charge) << "\n";
    }

    TH2F* hist = new TH2F(name.c_str(), name.c_str(), maxtick + 1, 0, maxtick + 1, maxchan + 1, 0, maxchan + 1);

    for (auto trace : (*traces)) {
        int tick = trace->tbin();
        int chan = trace->channel();
        auto& charge = trace->charge();
        for (float q : charge) {
            hist->Fill(tick + 0.5, chan + 0.5, q);
            ++tick;
        }
    }
    hist->SetXTitle("ticks");
    hist->SetYTitle("channels");
    hist->SetContour(50);
    hist->Draw("colz");
    if (mtick > 0) {
        hist->GetXaxis()->SetRangeUser(0, mtick);
    }
    if (mchan > 0) {
        hist->GetYaxis()->SetRangeUser(0, mchan);
    }

    pdf();
    return hist;
}

int main(int argc, char* argv[])
{
    PluginManager& pm = PluginManager::instance();
    auto aux_pi [[maybe_unused]] = pm.add("WireCellAux");
    assert(aux_pi);
    pm.add("WireCellGen");
    pm.add("WireCellRoot");

    auto idft = Factory::lookup_tn<IDFT>("FftwDFT");

    int nsamples = 50;
    double gain, shaping, tick;
    {
        IConfigurable::pointer mccfg = Factory::lookup<IConfigurable>("Misconfigure");
        Configuration cfg = mccfg->default_configuration();
        cfg["nsamples"] = nsamples;
        cfg["truncate"] = true;
        // use defaults locally
        gain = cfg["from"]["gain"].asDouble();
        shaping = cfg["from"]["shaping"].asDouble();
        tick = cfg["tick"].asInt();
        mccfg->configure(cfg);
    }
    IFrameFilter::pointer mc = Factory::lookup<IFrameFilter>("Misconfigure");

    Response::ColdElec ce(gain, shaping);

    auto resp = ce.generate(Binning(200, 0, 200 * tick));
    auto resp2 = ce.generate(Binning(400, 0, 400 * tick));
    auto resp3 = ce.generate(Binning(50, 0, 50 * tick));
    auto resp_spec = fwd_r2c(idft, resp);
    auto resp_spec2 = fwd_r2c(idft, resp2);
    auto resp_spec3 = fwd_r2c(idft, resp3);

    ITrace::vector q_traces;
    ITrace::vector out_traces;

    int qchannel = 0;
    int channel = 0;
    for (int nbins : {50, 100, 200, 400}) {
        Binning bins(nbins, 0, nbins * tick);

        auto q1 = blip(bins, 0.25 * nbins * tick);
        auto q2 = blip(bins, 0.50 * nbins * tick, 0.5 * units::us);
        auto q3 = blip(bins, 0.75 * nbins * tick);
        auto q4 = Waveform::realseq_t(nbins, 0);
        q4[nbins / 2] = maximum(q1);

        q_traces.push_back(std::make_shared<SimpleTrace>(qchannel++, 0, q1));
        q_traces.push_back(std::make_shared<SimpleTrace>(qchannel++, 0, q2));
        q_traces.push_back(std::make_shared<SimpleTrace>(qchannel++, 0, q3));
        q_traces.push_back(std::make_shared<SimpleTrace>(qchannel++, 0, q4));

        auto e1 = convolve(idft, q1, resp);
        auto e2 = convolve(idft, q2, resp);
        auto e3 = convolve(idft, q3, resp);
        auto e4 = convolve(idft, q4, resp);

        out_traces.push_back(std::make_shared<SimpleTrace>(channel++, 0, e1));
        out_traces.push_back(std::make_shared<SimpleTrace>(channel++, 0, e2));
        out_traces.push_back(std::make_shared<SimpleTrace>(channel++, 0, e3));
        out_traces.push_back(std::make_shared<SimpleTrace>(channel++, 0, e4));
    }

    auto qorig = std::make_shared<SimpleFrame>(0, 0.0, q_traces, tick);
    auto orig = std::make_shared<SimpleFrame>(0, 0.0, out_traces, tick);
    IFrame::pointer mced, dummy;
    {
        bool ok = (*mc)(orig, mced);
        if (!ok) {
            return -1;
        }
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    MultiPdf pdf(argv[0]);

    plot_wave(pdf, resp, "resp", "Electronics response waveform - 200 bins");
    plot_wave(pdf, resp2, "resp2", "Electronics response waveform - 400 bins");
    plot_wave(pdf, resp3, "resp3", "Electronics response waveform - 50 bins");
    plot_wave(pdf, Waveform::real(resp_spec), "respect", "Electronics response spectrum - 200 bins");
    plot_wave(pdf, Waveform::real(resp_spec2), "respect2", "Electronics response spectrum - 400 bins");
    plot_wave(pdf, Waveform::real(resp_spec3), "respect3", "Electronics response spectrum - 50 bins");

    plot_frame(pdf, qorig, "ChargeZoomed", 50., 5.);
    plot_frame(pdf, orig, "ShapedZoomed", 50., 5.);
    plot_frame(pdf, mced, "MisconfiguredZoomed", 50, 5);

    plot_frame(pdf, qorig, "Charge");
    auto f1 = plot_frame(pdf, orig, "Shaped");
    auto f2 = plot_frame(pdf, mced, "Misconfigured");
    plot_ratio(pdf, f1, f2);

    return 0;
}
