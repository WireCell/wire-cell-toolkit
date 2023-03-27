#include "WireCellAux/TensorDM.h"
#include "WireCellAux/SimpleTensorSet.h"
#include <regex>

using namespace WireCell;

ITensorSet::pointer
WireCell::Aux::TensorDM::as_tensorset(const ITensor::vector& tens,
                                      int ident,
                                      const Configuration& tsetmd)
{
    auto sv = std::make_shared<ITensor::vector>(tens.begin(), tens.end());
    return std::make_shared<SimpleTensorSet>(ident, tsetmd, sv);
}


WireCell::Aux::TensorDM::located_t
WireCell::Aux::TensorDM::index_datapaths(const ITensor::vector& tens)
{
    WireCell::Aux::TensorDM::located_t located;
    for (const auto& iten : tens) {
        auto md = iten->metadata();
        auto jdp = md["datapath"];
        if (!jdp.isString()) {
            continue;
        }
        auto dp = jdp.asString();
        located[dp] = iten;
    }
    return located;
}

ITensor::pointer
WireCell::Aux::TensorDM::first_of(const ITensor::vector& tens,
                                  const std::string& datatype)
{
    for (const auto& iten : tens) {
        auto md = iten->metadata();
        auto jdt = md["datatype"];
        if (!jdt.isString()) {
            continue;
        }
        auto dt = jdt.asString();
        if (dt == datatype) {
            return iten;
        }
    }
    return nullptr;
}
ITensor::vector
WireCell::Aux::TensorDM::match_at(const ITensor::vector& tens,
                                  const std::string& datapath)
{
    std::regex re(datapath);
    ITensor::vector ret;
    for (const auto& iten : tens) {
        auto md = iten->metadata();
        auto jdp = md["datapath"];
        if (!jdp.isString()) {
            continue;
        }
        auto dp = jdp.asString();
        if (! std::regex_match(dp, re)) continue;
        ret.push_back(iten);
    }
    return ret;
}

