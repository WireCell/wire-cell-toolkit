/** Produce a schema of all registered components.

    Eg

    wire-cell -a ConfigSchema -p WireCellApps -p WireCellGen

 */

#ifndef WIRECELLAPPS_CONFIGSCHEMA
#define WIRECELLAPPS_CONFIGSCHEMA

#include "WireCellIface/IApplication.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Configuration.h"

namespace WireCellApps {

    class ConfigSchema : public WireCell::IApplication, public WireCell::IConfigurable {
        WireCell::Configuration m_cfg;

       public:
        ConfigSchema();
        virtual ~ConfigSchema();

        virtual void execute();

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
    };

}  // namespace WireCellApps

#endif
