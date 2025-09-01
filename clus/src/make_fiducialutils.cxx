/// This is an ensemble visitor that constructs a FiducialUtils and places it in
/// the "live" grouping.

#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"

#include <memory> // for make_shared


class MakeFiducialUtils;
WIRECELL_FACTORY(MakeFiducialUtils, MakeFiducialUtils,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)
using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
class MakeFiducialUtils : public IConfigurable, public Clus::IEnsembleVisitor
                        , private Clus::NeedDV, private Clus::NeedFiducial, private Clus::NeedPCTS {
public:
    MakeFiducialUtils() {};
    virtual ~MakeFiducialUtils() {};

    void configure(const WireCell::Configuration& cfg) {
        NeedDV::configure(cfg);
        NeedFiducial::configure(cfg);
        NeedPCTS::configure(cfg);
        m_live_name = get(cfg, "live", m_live_name);
        m_dead_name = get(cfg, "dead", m_dead_name);
        m_targ_name = get(cfg, "target", m_targ_name);
        

    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["live"] = m_live_name;
        cfg["dead"] = m_dead_name;
        cfg["target"] = m_targ_name;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        using spdlog::warn;

        auto* live_grouping = ensemble.with_name(m_live_name).at(0);
        auto* dead_grouping = ensemble.with_name(m_dead_name).at(0);
        auto* targ_grouping = ensemble.with_name(m_targ_name).at(0);
        
        int nerrors = 0;
        if (!live_grouping) {
            warn("no live grouping in ensemble");
            ++nerrors;
        }
        if (!dead_grouping) {
            warn("no dead grouping in ensemble");
            ++nerrors;
        }
        if (!targ_grouping) {
            warn("no target grouping in ensemble");
            ++nerrors;
        }
        if (nerrors) {
            warn("ensemble is not well formed, no FiducialUtils made, this may cause downstream problems.");
            return;
        }
        
        std::cout << "Add FiducialUtils to grouping" << std::endl;
        auto fu = std::make_shared<FiducialUtils>(FiducialUtils::StaticData{m_dv, m_fiducial, m_pcts});

        // Feed the dynamic data (live and dead groupings)
        fu->feed_dynamic(FiducialUtils::DynamicData{*live_grouping, *dead_grouping});
        
        // validation ...
        // {
        // // fu->inside_fiducial_volume(Point(250*units::cm, 120*units::cm, 10*units::cm));
        // // fu->inside_fiducial_volume(Point(250*units::cm, -80*units::cm, 30*units::cm));
        // // fu->inside_fiducial_volume(Point(250*units::cm, 120*units::cm, -10*units::cm));
        // // fu->inside_fiducial_volume(Point(25*units::cm, 60*units::cm, 10*units::cm));
        // // fu->inside_fiducial_volume(Point(250*units::cm, 120*units::cm, 10*units::cm));
        // // fu->inside_fiducial_volume(Point(250*units::cm, -80*units::cm, 1030*units::cm));
        // // fu->inside_fiducial_volume(Point(250*units::cm, -80*units::cm, 960*units::cm));
        // // fu->inside_fiducial_volume(Point(25*units::cm, 60*units::cm, 10*units::cm));

        //     std::cout << "FiducialUtils check: " << fu->inside_dead_region(Point(0,0,2),0,0,2) << " " << fu->inside_dead_region(Point(41.8*units::cm, 26.5*units::cm,707.0*units::cm),0,0,2) << " " << fu->inside_dead_region(Point(41.8*units::cm, 16.5*units::cm,707.0*units::cm),0,0,2) << std::endl; 


        // }

        live_grouping->set_fiducialutils(fu);
    }

private:

    // Names of grouping for "live" and "dead" and "target" (the grouping to
    // receive the FiducialUtils.
    std::string m_live_name{"live"}, m_dead_name{"dead"}, m_targ_name{"live"};
};
