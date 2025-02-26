#ifndef WIRECELLAPPS_CONFIGDUMPER
#define WIRECELLAPPS_CONFIGDUMPER

#include "WireCellIface/IApplication.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Configuration.h"

namespace WireCellApps {

    class AnodeDumper : public WireCell::IApplication, public WireCell::IConfigurable {
        WireCell::Configuration m_cfg;

       public:
        AnodeDumper();
        virtual ~AnodeDumper();

        virtual void execute();

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
    };

}  // namespace WireCellApps

#endif
