#ifndef WIRECELLAPPS_SCHEMADUMPER
#define WIRECELLAPPS_SCHEMADUMPER

#include "WireCellIface/IApplication.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Configuration.h"

namespace WireCellApps {

    class SchemaDumper : public WireCell::IApplication, public WireCell::IConfigurable {
        WireCell::Configuration m_cfg;

       public:
        SchemaDumper();
        virtual ~SchemaDumper();

        virtual void execute();

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
    };

}  // namespace WireCellApps

#endif
