/** A base class providing auto-configuration.

    This works in conjunction with a config struct which is produced
    via shema codegen.

    To use

    1) Inherit your MyConfigurable class from
    Configurable<Cfg::MyConfigurable> instead of directly
    IConfigurable.

    2) If migrating from old, remove configure() and
    default_configure() methods.

    3) If you configure() method held non-configuration related
    initialization, it may move into configured().

    4) Use the protected m_cfg object to access the configuration
    struct including any "service" type components "used" by yours.

 */


#ifndef WIRECELLAUX_CONFIGURABLE
#define WIRECELLAUX_CONFIGURABLE

#include "WireCellUtil/nljs2jcpp.hpp" // remove when ditch JsonCPP

#include "WireCellIface/IConfigurable.h"

namespace WireCell::Aux {
    
    template<typename CfgStruct>
    class Configurable : public IConfigurable {
      public:
        using cfg_type = CfgStruct;

        // Subclass may implement in order to get post-config entry.
        virtual void configured() { }

        // Subclass need not but MAY overide IF it also calls.
        virtual void configure(const WireCell::Configuration& jcfg)
        {
            nljs_t nljs = jcfg;
            m_cfg = nljs.get<cfg_type>();
            configured();
        }
        
        // Subclass need not but MAY overide IF it also calls.
        virtual WireCell::Configuration default_configuration() const
        {
            const nljs_t nljs = m_cfg;
            return nljs.get<Json::Value>();
        }

      protected:
        // Subclass should use m_cfg directly.
        cfg_type m_cfg;
    };

}


#endif
