
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"

#include "WireCellUtil/Exceptions.h"

#include <vector>
#include <string>

class CompositeFiducial;

WIRECELL_FACTORY(CompositeFiducial, CompositeFiducial,
                 WireCell::IFiducial, WireCell::IConfigurable)

using namespace WireCell;


class CompositeFiducial : public WireCell::IFiducial, public WireCell::IConfigurable {
public:
    CompositeFiducial();
    virtual ~CompositeFiducial();

    // IConfigurable
    virtual void configure(const WireCell::Configuration& cfg);
    virtual WireCell::Configuration default_configuration() const;

    // IFiducial
    virtual bool contained(const WireCell::Point& point) const;

private:
    // Vector of child IFiducial components
    std::vector<WireCell::IFiducial::pointer> m_fiducials;
    
    // Logic for combining results: "and", "or", "nand", "nor"
    std::string m_logic{"and"};
    
    // Helper method to apply logic
    bool apply_logic(const std::vector<bool>& results) const;
};




CompositeFiducial::CompositeFiducial() {}

CompositeFiducial::~CompositeFiducial() {}

Configuration CompositeFiducial::default_configuration() const {
    Configuration cfg;
    
    // Logic for combining child fiducial results
    cfg["logic"] = m_logic;
    
    // Array of child IFiducial component names to combine
    cfg["fiducials"] = Json::arrayValue;
    
    return cfg;
}

void CompositeFiducial::configure(const Configuration& cfg) {
    m_logic = get(cfg, "logic", m_logic);
    
    // Validate logic
    if (m_logic != "and" && m_logic != "or" && m_logic != "nand" && m_logic != "nor") {
        raise<ValueError>("CompositeFiducial: invalid logic '%s', must be 'and', 'or', 'nand', or 'nor'", m_logic);
    }
    
    // Clear existing fiducials
    m_fiducials.clear();
    
    // Load child fiducial components
    const auto& fiducial_names = cfg["fiducials"];
    if (fiducial_names.empty()) {
        raise<ValueError>("CompositeFiducial: 'fiducials' array cannot be empty");
    }
    
    // std::cout << "CompositeFiducial: Loading " << fiducial_names.size() << " child fiducials with logic '" << m_logic << "'" << std::endl;
    
    for (const auto& name : fiducial_names) {
        std::string fiducial_tn = name.asString();
        // std::cout << "CompositeFiducial: Looking for fiducial '" << fiducial_tn << "'" << std::endl;
        auto fiducial = Factory::find_tn<IFiducial>(fiducial_tn);
        if (!fiducial) {
            raise<ValueError>("CompositeFiducial: failed to find IFiducial component '%s'", fiducial_tn);
        }
        // std::cout << "CompositeFiducial: Successfully loaded fiducial '" << fiducial_tn << "'" << std::endl;
        m_fiducials.push_back(fiducial);
    }
}

bool CompositeFiducial::contained(const Point& point) const {
    if (m_fiducials.empty()) {
        return true;  // No restrictions if no child fiducials
    }
    
    // Evaluate all child fiducials
    std::vector<bool> results;
    results.reserve(m_fiducials.size());
    
    for (size_t i = 0; i < m_fiducials.size(); ++i) {
        const auto& fiducial = m_fiducials[i];
        bool result = fiducial->contained(point);
        results.push_back(result);
        // std::cout << "CompositeFiducial: point " << point << " contained by " 
                //   << (result ? "true" : "false") << " in fiducial[" << i << "]" << std::endl;
    }
    
    // Apply combination logic
    return apply_logic(results);
}

bool CompositeFiducial::apply_logic(const std::vector<bool>& results) const {
    if (results.empty()) return true;
    
    if (m_logic == "and") {
        // All must be true
        for (bool result : results) {
            if (!result) return false;
        }
        return true;
    }
    else if (m_logic == "or") {
        // At least one must be true
        for (bool result : results) {
            if (result) return true;
        }
        return false;
    }
    else if (m_logic == "nand") {
        // NOT(all true) = at least one false
        for (bool result : results) {
            if (!result) return true;
        }
        return false;
    }
    else if (m_logic == "nor") {
        // NOT(any true) = all false
        for (bool result : results) {
            if (result) return false;
        }
        return true;
    }
    
    return true;  // Default fallback
}