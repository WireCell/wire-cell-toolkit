#include "WireCellTrig/TASetFanin.h"

#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>

WIRECELL_FACTORY(TASetFanin, WireCell::Trig::TASetFanin, WireCell::INamed, WireCell::ITriggerActivitySetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Trig;

TASetFanin::TASetFanin(size_t multiplicity)
  : Aux::Logger("TASetFanin", "trig")
  , m_multiplicity(multiplicity)
{
}

TASetFanin::~TASetFanin() {}

WireCell::Configuration TASetFanin::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = (int) m_multiplicity;
    return cfg;
}

void TASetFanin::configure(const WireCell::Configuration& cfg)
{
    const int mult = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (mult <= 0) {
        THROW(ValueError() << errmsg{"TASetFanin: multiplicity must be positive"});
    }
    m_multiplicity = mult;
}

std::vector<std::string> TASetFanin::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

bool TASetFanin::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    if (invec.size() != m_multiplicity) {
        THROW(ValueError() << errmsg{"TASetFanin: got " + std::to_string(invec.size()) + " inputs, expected " +
                                     std::to_string(m_multiplicity)});
    }
    const size_t neos = std::count(invec.begin(), invec.end(), nullptr);
    if (neos) {
        if (neos != invec.size()) {
            log->warn("end of stream on {} of {} ports at call {}", neos, invec.size(), m_count);
        }
        return true;
    }

    const auto& first = invec.front();
    const double begin = first->time_begin(), end = first->time_end();
    const double tol = 1 * units::ns;
    int detector = first->detector();
    ITriggerActivitySet::activities_t tas;
    for (const auto& one : invec) {
        if (std::abs(one->time_begin() - begin) > tol || std::abs(one->time_end() - end) > tol) {
            THROW(ValueError() << errmsg{"TASetFanin: input coverages differ at call " + std::to_string(m_count)});
        }
        if (one->detector() != detector) detector = -1;
        const auto& more = one->activities();
        tas.insert(tas.end(), more.begin(), more.end());
    }
    std::stable_sort(tas.begin(), tas.end(),
                     [](const auto& a, const auto& b) { return a.time_start < b.time_start; });
    out = std::make_shared<Aux::SimpleTriggerActivitySet>(m_count++, detector, begin, end, std::move(tas));
    return true;
}
