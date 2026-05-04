#include "WireCellRoot/HistFrameSink.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Waveform.h"

#include "TFile.h"
#include "TH2F.h"

#include <iostream>
#include <algorithm>
#include <unordered_map>

WIRECELL_FACTORY(HistFrameSink, WireCell::Root::HistFrameSink, WireCell::IFrameSink, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Root::HistFrameSink::HistFrameSink()
  : m_filepat("histframe-%02d.root")
  , m_anode_tn("AnodePlane")
  , m_units(1.0)
{
}

Root::HistFrameSink::~HistFrameSink() {}

WireCell::Configuration Root::HistFrameSink::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filepat;
    cfg["anode"] = m_anode_tn;
    cfg["units"] = units::mV;
    return cfg;
}

void Root::HistFrameSink::configure(const WireCell::Configuration& cfg)
{
    m_filepat = get<std::string>(cfg, "filename", m_filepat);
    m_units = get<double>(cfg, "units", m_units);
    m_anode_tn = get<std::string>(cfg, "anode", m_anode_tn);
    m_anode = Factory::lookup_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        cerr << "Root::HistFrameSink: failed to get anode: \"" << m_anode_tn << "\"\n";
        return;
    }
    cerr << "Root::HistFrameSink: configured with: "
         << "file:" << m_filepat << ", "
         << "units:" << m_units << ", "
         << "anode:" << m_anode_tn << endl;
}

bool Root::HistFrameSink::operator()(const IFrame::pointer& frame)
{
    if (!frame) {
        cerr << "Root::HistFrameSink: no frame\n";
        return true;
    }

    std::string fname = Form(m_filepat.c_str(), frame->ident());
    TFile* file = TFile::Open(fname.c_str(), "recreate");

    typedef std::tuple<ITrace::vector, std::vector<int>, std::vector<int> > tct_tuple;
    std::unordered_map<int, tct_tuple> perplane;

    // collate traces into per plane and also calculate bounds
    ITrace::shared_vector traces = frame->traces();
    for (ITrace::pointer trace : *traces) {
        int ch = trace->channel();
        auto wpid = m_anode->resolve(ch);
        int wpident = wpid.ident();
        double tmin = trace->tbin();
        double tlen = trace->charge().size();

        auto& tct = perplane[wpident];
        get<0>(tct).push_back(trace);
        get<1>(tct).push_back(ch);
        get<2>(tct).push_back(tmin);
        get<2>(tct).push_back(tmin + tlen);

        if (wpident < 0) {
            cerr << "Channel " << ch << " has illegal wire plane ident: " << wpid << endl;
        }
    }

    const double t0 = frame->time();
    const double tick = frame->tick();

    for (auto& thisplane : perplane) {
        int wpident = thisplane.first;
        auto& tct = thisplane.second;
        auto& traces = get<0>(tct);
        auto& chans = get<1>(tct);
        auto chmm = std::minmax_element(chans.begin(), chans.end());
        auto& tbins = get<2>(tct);
        auto tbmm = std::minmax_element(tbins.begin(), tbins.end());

        const double tmin = t0 + tick * (*tbmm.first);
        const double tmax = t0 + tick * (*tbmm.second);
        const int ntbins = (*tbmm.second) - (*tbmm.first);

        const int chmin = round(*chmm.first);
        const int chmax = round(*chmm.second + 1);
        const int nchbins = chmax - chmin;

        TH2F* hist = new TH2F(Form("plane%d", wpident), Form("Plane %d", wpident), ntbins, tmin / units::us,
                              tmax / units::us, nchbins, chmin, chmax);
        hist->SetDirectory(file);  // needed?
        hist->SetXTitle("time [us]");
        hist->SetYTitle("channel");

        const int tbin_offset = *tbmm.first;  // min tbin across all traces in this plane
        double qtot = 0;
        int nbins_tot = 0;
        for (auto& trace : traces) {
            int ch = trace->channel();
            int tbin = trace->tbin();
            auto& charge = trace->charge();
            int nbins = charge.size();
            nbins_tot += nbins;
            int ybin = ch - chmin + 1;  // ROOT bin (1-based)
            for (int ibin = 0; ibin < nbins; ++ibin) {
                int xbin = (tbin + ibin) - tbin_offset + 1;  // ROOT bin (1-based)
                double val = charge[ibin] / m_units;
                hist->SetBinContent(xbin, ybin, val + hist->GetBinContent(xbin, ybin));
                qtot += charge[ibin];
            }
        }

        cerr << wpident << " ntraces:" << traces.size() << " "
             << " nsamples:" << nbins_tot << " "
             << " qtot:" << qtot / m_units << " "
             << " qunit:" << m_units << " "
             << " integ:" << hist->Integral() << " min:" << hist->GetMinimum() << " max:" << hist->GetMaximum()
             << " chan:" << nchbins << "[" << chmin << "," << chmax << "] "
             << " time:" << ntbins << "[" << tmin / units::us << "," << tmax / units::us << "]us\n";

        hist->Write();
    }
    file->Close();
    delete file;
    file = nullptr;
    return true;
}
