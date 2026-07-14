#include "WireCellRoot/RootfileCreation.h"

#include "TFile.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(RootfileCreation_depos, WireCell::Root::RootfileCreation_depos, WireCell::IDepoFilter,
                 WireCell::IConfigurable)

WIRECELL_FACTORY(RootfileCreation_frames, WireCell::Root::RootfileCreation_frames, WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

Root::RootfileCreation_depos::RootfileCreation_depos() {}

Root::RootfileCreation_depos::~RootfileCreation_depos() {}

void Root::RootfileCreation_depos::configure(const WireCell::Configuration& cfg)
{
    m_cfg = cfg;
    create_file();
}

WireCell::Configuration Root::RootfileCreation_depos::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = "";
    cfg["root_file_mode"] = "RECREATE";
    return cfg;
}

bool Root::RootfileCreation_depos::operator()(const WireCell::IDepo::pointer& indepo, WireCell::IDepo::pointer& outdepo)
{
    outdepo = indepo;
    if (!indepo) {
        // eos
        std::cerr << "RootfileCreation_depos: EOS\n";
        return true;
    }
    return true;
}

void Root::RootfileCreation_depos::create_file()
{
    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = m_cfg["root_file_mode"].asString();
    TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());
    output_tf->Close("R");
    delete output_tf;
    output_tf = nullptr;
}

Root::RootfileCreation_frames::RootfileCreation_frames() {}

Root::RootfileCreation_frames::~RootfileCreation_frames() {}

void Root::RootfileCreation_frames::configure(const WireCell::Configuration& cfg)
{
    m_cfg = cfg;
    create_file();
}

WireCell::Configuration Root::RootfileCreation_frames::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = "";
    cfg["root_file_mode"] = "RECREATE";
    return cfg;
}

void Root::RootfileCreation_frames::create_file()
{
    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = m_cfg["root_file_mode"].asString();

    //  std::cout << "Xin: " << ofname << std::endl;

    TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());
    output_tf->Close("R");
    delete output_tf;
    output_tf = nullptr;
}

bool Root::RootfileCreation_frames::operator()(const WireCell::IFrame::pointer& in, WireCell::IFrame::pointer& out)
{
    out = in;
    if (!in) {
        // eos
        std::cerr << "RootfileCreation_frames: EOS\n";
        return true;
    }
    return true;
}
