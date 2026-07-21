#pragma once

/**
   This API provides helpers to deal with aggregating configurable classes
   either through inheritance or member variables.

   You class may be:

   - itself configurable with JSON or C++ struct (via Hana)
   - inherit from zero or more configurable classes.
   - have zero or more configurable member variables.

   With the following constraints:

   - No non-virtual diamond inheritance pattern.

   For illustration, assume your class MyClass is to be configurable, and
   inherits from OtherConfigurable and has ExtraConfigurable as a member
   variable.

   It is best that you let Hana handle the serializing the Configurable JSON
   object to/from a configuration struct for your component like this:

   @code{cpp}
   struct MyConfig {};
   struct MyClass : HanaConfigurable<MyClass, MyConfig, OtherConfigurable> {
       ExtraConfigurable extra;
 
       virtual WireCell::Configuration default_configuration_extra() const {
           return extra.default_configuration();
       }
       virtual void configure_extra(const WireCell::Configuration& jconfig) {
           extra.configure(jconfig);
       }

       virtual void configured() {
           // use "m_config" variable of type MyConfig
       }
   };
   @endcode

   If your class has no configurable data members, there is no need to supply
   the `*_extra()` methods.  If it has more than 1, see the `*_each()` functions
   in the next example.

   Not recommended but if you really want to directly handle JSON serialization
   you can do so by inheriting from AutoConfigurable like:

    @code{cpp}
    struct MyClass : AutoConfigurable<MyClass, OtherConfigurable> : public virtual IConfigurable {
        ExtraConfigurable extra;
        YetMoreConfigurable yetmore;

        MyClass() {
            // Can provide defaults into the protected variable:
            m_json_config["myparam"] = 42;
        }
        
        virtual WireCell::Configuration default_configuration_json() const {
             Configuration cfg = default_configuration_each({&extra,&yetmore});
             // You could also provide defaults here in "cfg" but they are not
             // stored into m_json_config.
             return cfg;
        }

        virtual void configured() {
            // Here we get notification that our m_json_config has been set.

            // Dispatch to members.  Dispatch to parent class is automatic.
            configure_each({&extra,&yetmore}, m_json_config);

            // Interpret the jconfig for MyClass's own parts.
        }
    };
    @endcode

    Here, we also show the use of default_configuration_each() and
    configure_each() which simply iterates over its arguments to call the
    IConfigurable methods.


    If your class merely inherits from configurables and wants to provide a
    unified IConfigurable interface, it may inherit from AutoConfigurable and
    simply not provide any configuration related methods.

 */


#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellIface/IConfigurable.h"


// This could move into Aux so I put it in a more generic namespace
namespace WireCell {

    /// Helper to configure a parameter pack of different types.
    template<typename Type>
    void configure_types(const WireCell::Configuration& config, Type& obj) {
        obj.configure(config);
    }
    template<typename Type, typename... Types>
    void configure_types(const WireCell::Configuration& config, Type& obj, Types&... rest) {
        obj.configure(config);
        configure_types(config, rest...);
    }

    /// Helper to aggregate default config over many.
    template<typename Type>
    WireCell::Configuration default_configuration_types(const Type& obj) {
        return obj.default_configuration();
    }
    template<typename Type, typename... Types>
    WireCell::Configuration default_configuration_types(const Type& obj, const Types&... rest)  {
        WireCell::Configuration cfg = obj.default_configuration();
        WireCell::Configuration cfg2 = default_configuration_types(rest...);
        return update(cfg, cfg2);
    }

    template<typename ThisType, typename... Bases>
    void configure_bases(ThisType* this_obj, const Configuration& cfg)  {
        (void)(((static_cast<Bases*>(static_cast<ThisType*>(this_obj)))->Bases::configure(cfg)), ...);
    }

    template<typename ThisType, typename... Bases>
    WireCell::Configuration default_configuration_bases(const ThisType* this_obj)  {
        WireCell::Configuration aggregated_cfg;

        (void)(((
                    aggregated_cfg = update(
                        aggregated_cfg,
                        (static_cast<const Bases*>(static_cast<const ThisType*>(this_obj)))->Bases::default_configuration()
                        )
                    )), ...);

        return aggregated_cfg;

    }
    

    // Helper to aggregate default configuration across many.
    inline
    WireCell::Configuration default_configuration_each(std::initializer_list<const IConfigurable*> objs)
    {
        WireCell::Configuration combined_config;

        // We iterate over the list. T will be deduced based on what is passed (e.g., A&, B&, C&).
        for (const auto& obj : objs) {
            auto current_config = obj->default_configuration();
            update(combined_config, current_config);
        }
        return combined_config;
    }

    // Helper to dispatch configuration to many.
    inline
    void configure_each(std::initializer_list<IConfigurable*> objs, const WireCell::Configuration& config)
    {
        for (auto& obj : objs) {
            obj->configure(config);
        }
    }




    // Helper to detect a method.
    template <typename T, typename = void>
    struct has_member_configured : std::false_type {};
    template <typename T>
    struct has_member_configured<T, std::void_t<decltype(&T::configured)>> : std::true_type {};

    /// Inherit like
    /// struct AC : AutoConfigurable<AC, Base1, Base2>
    template <typename Derived, typename... Bases>
    struct AutoConfigurable : public Bases..., virtual public IConfigurable {

        template <typename Parent>
            void call_configured() {
            if constexpr (has_member_configured<Parent>::value) {
                static_cast<Parent*>(this)->Parent::configured();
            }
        }

        // Derived may override to get the configuration in m_json_config.
        virtual void configured() { }

        // Derived may override to provide its default.
        virtual WireCell::Configuration default_configuration_json() const {
            return m_json_config;
        }

        // Dispatch configuration to base classes.
        void configure_bases() {
            (void)(((static_cast<Bases*>(static_cast<Derived*>(this)))->Bases::configure(m_json_config)), ...);
        }

        void configured_bases() {
            (call_configured<Bases>(), ...);
        }

        // Collect default configuration across base classes.
        WireCell::Configuration default_configuration_bases() const {
            WireCell::Configuration aggregated_cfg;

            (void)(((
                        aggregated_cfg = update(
                            aggregated_cfg,
                            (static_cast<const Bases*>(static_cast<const Derived*>(this)))->Bases::default_configuration()
                            )
                        )), ...);

            return aggregated_cfg;
        }

        // Dispatch to base classes and to self.
        virtual void configure(const WireCell::Configuration& jconfig) {
            m_json_config = jconfig;
            configure_bases();
            configured_bases();
        }

        // Aggregate across mixins and self.
        virtual WireCell::Configuration default_configuration() const {
            auto cfg = default_configuration_bases();
            auto cfg2 = default_configuration_json();
            update(cfg, cfg2);
            return cfg;
        }

    protected:
        Configuration m_json_config;

    };
    // use, eg
    // struct Threshold : public AutoConfigurable<Threshold, ContextBase, Logger> { /*...*/ };


    template<typename Derived, typename Config, typename... Bases>
    struct HanaConfigurable : public Bases..., virtual public IConfigurable
    {

        template <typename Parent>
            void call_configured() {
            if constexpr (has_member_configured<Parent>::value) {
                static_cast<Parent*>(this)->Parent::configured();
            }
        }

        void configure_bases(const WireCell::Configuration& config) {
            // C++17 fold expression over the comma operator
            // This expands to:
            // (static_cast<Bases*>(static_cast<Derived*>(this))->configure(config), ...);
            // The 'void()' ensures it's a statement.
            (void)(((static_cast<Bases*>(static_cast<Derived*>(this)))->Bases::configure(config)), ...);
        }

        void configured_bases() {
            (call_configured<Bases>(), ...);
        }

        virtual WireCell::Configuration default_configuration_extra() const {
            Configuration cfg;
            return cfg;
        }

        WireCell::Configuration default_configuration_bases() const {
            WireCell::Configuration aggregated_cfg;

            // C++17 fold expression over the comma operator.
            // The lambda captures 'aggregated_cfg' by reference and updates it.
            (void)(((
                        aggregated_cfg = update(
                            aggregated_cfg,
                            (static_cast<const Bases*>(static_cast<const Derived*>(this)))->Bases::default_configuration()
                            )
                        )), ...);

            return aggregated_cfg;
        }

        // Derived may override to get notification that m_hana_config has been set.
        virtual void configured() { }

        // Get notification of the JSON configuration in order to configure any
        // parts not covered by Hana.  For example, if your class holds other
        // configurables as member variables.
        virtual void configure_extra(const WireCell::Configuration& jconfig) {
            return;
        }

        // Dispatch to base mixins and to self.
        virtual void configure(const WireCell::Configuration& jconfig) {
            configure_bases(jconfig);
            configure_extra(jconfig);
            HanaJsonCPP::from_json(m_hana_config, jconfig);
            configured_bases();
            configured();
        }

        // Aggregate across mixins and self.
        virtual WireCell::Configuration default_configuration() const {
            auto cfg = default_configuration_bases();
            auto cfg2 = HanaJsonCPP::to_json(m_hana_config);
            update(cfg, cfg2);
            cfg2 = default_configuration_extra();
            update(cfg, cfg2);            
            return cfg;
        }

        Config& hana_config() { return m_hana_config; }
        const Config& hana_config() const { return m_hana_config; }


    protected:

        Config m_hana_config;
    };


}
