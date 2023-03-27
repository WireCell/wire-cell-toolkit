#include "WireCellAux/FrameTensor.h"
#include "WireCellAux/TensorDM.h"

#include "WireCellUtil/String.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(FrameTensor, WireCell::Aux::FrameTensor,
                 WireCell::IFrameTensorSet,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;

FrameTensor::FrameTensor()
{
}
FrameTensor::~FrameTensor()
{
}


bool FrameTensor::operator()(const input_pointer& frame, output_pointer& its)
{
    its = nullptr;
    if (! frame) {
        // eos
        return true;
    }

    const int ident = frame->ident();
    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }

    auto tens = as_tensors(frame, datapath, m_mode, m_digitize, m_transform);
    its = as_tensorset(tens, ident);
    return true;
}
        


WireCell::Configuration FrameTensor::default_configuration() const
{
    std::map<FrameTensorMode, std::string> modes = {
        {FrameTensorMode::unified, "unified"},
        {FrameTensorMode::tagged, "tagged" },
        {FrameTensorMode::sparse, "sparse" }};
    Configuration cfg;
    cfg["mode"] = modes[m_mode];
    cfg["digitize"] = m_digitize;
    cfg["baseline"] = m_baseline;
    cfg["offset"] = m_offset;
    cfg["scale"] = m_scale;
    cfg["datapath"] = m_datapath;

    return cfg;
}
void FrameTensor::configure(const WireCell::Configuration& cfg)
{
    m_datapath = get(cfg, "datapath", m_datapath);
    m_baseline = get(cfg, "baseline", m_baseline);
    m_offset = get(cfg, "offset", m_offset);
    m_scale = get(cfg, "scale", m_scale);
    m_transform = [&](float x) -> float {
        return (x + m_baseline) * m_scale + m_offset;
    };

    std::map<std::string, FrameTensorMode> modes = {
        {"unified",  FrameTensorMode::unified},
        {"tagged", FrameTensorMode::tagged },
        {"sparse", FrameTensorMode::sparse }};
    auto mode = get<std::string>(cfg, "mode", "tagged");
    m_mode = modes[mode];
    m_digitize = get(cfg, "digitize", m_digitize);
}

        
