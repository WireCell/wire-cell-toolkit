#include "WireCellAux/TensorFrame.h"
#include "WireCellAux/TensorDM.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(TensorFrame, WireCell::Aux::TensorFrame,
                 WireCell::ITensorSetFrame)

using namespace WireCell;
using namespace WireCell::Aux;

TensorFrame::TensorFrame()
{
}
TensorFrame::~TensorFrame()
{
}


const std::string default_dpre = "frames/[[:digit:]]+$";

WireCell::Configuration TensorFrame::default_configuration() const
{
    Configuration cfg;
    cfg["dpre"] = default_dpre;

    cfg["baseline"] = m_baseline;
    cfg["offset"] = m_offset;
    cfg["scale"] = m_scale;

    return cfg;
}
void TensorFrame::configure(const WireCell::Configuration& cfg)
{
    auto dpre = get(cfg, "dpre", default_dpre);
    m_dpre = std::regex(dpre);
    
    m_baseline = get(cfg, "baseline", m_baseline);
    m_offset = get(cfg, "offset", m_offset);
    m_scale = get(cfg, "scale", m_scale);

    m_transform = [&](float x) -> float {
        return (x - m_offset) / m_scale - m_baseline;
    };
}


bool TensorFrame::operator()(const input_pointer& its, output_pointer& frame)
{
    frame = nullptr;
    if (! its) {
        // eos
        return true;
    }

    const ITensor::vector& tens = *(its->tensors().get());
    for (const auto& ten : tens) {
        std::string dp = ten->metadata()["datapath"].asString();
        if (std::regex_match(dp, m_dpre)) {
            frame = TensorDM::as_frame(tens, dp, m_transform);
            break;
        }
    }
    return true;
}
