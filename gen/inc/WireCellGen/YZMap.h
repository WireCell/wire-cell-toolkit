#ifndef WIRECELL_GEN_YZMAP
#define WIRECELL_GEN_YZMAP

#include "WireCellIface/IYZMap.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace WireCell::Gen {

    class YZMap : public Aux::Logger, public IYZMap, public IConfigurable {
    public:
        YZMap();
        virtual ~YZMap();

        virtual double value(const std::string& anode_name,
                             int plane,
                             double y, double z) const;

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    private:
        // Full 4-level map: anode_name -> plane -> [binz][biny]
        std::unordered_map<std::string,
            std::unordered_map<int,
                std::vector<std::vector<double>>>> m_map;

        double m_bin_width{0};
        double m_bin_height{0};
        double m_yoffset{0};
        double m_zoffset{0};
        int m_nbinsy{0};
        int m_nbinsz{0};
    };

}

#endif
