#include <vector>
#ifndef WIRECELL_NAMEDFACTORY_H
#define WIRECELL_NAMEDFACTORY_H

#include "WireCellUtil/IFactory.h"
#include "WireCellUtil/Singleton.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Type.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include <unordered_map>
#include <functional>
#include <mutex>

#include <iostream>  // fixme: remove
#include <exception>
#include <string>
#include <set>

namespace WireCell {

    struct FactoryException : virtual public Exception {
    };

    /** Global registry that tracks all interface-specific factory registries.
     * This allows introspection of all registered factories without
     * compile-time knowledge of interface types.
     */
    class GlobalFactoryRegistry {
       public:
        /// Information about a registered interface type
        struct InterfaceInfo {
            std::string interface_name;  ///< Demangled interface type name

            /// Get all known type names for this interface
            std::function<std::vector<std::string>()> get_known_types;

            /// Instantiate an object by classname, returns Interface::pointer or nullptr
            std::function<Interface::pointer(const std::string&, bool)> instantiate;
        };

        /// Register an interface type with the global registry
        void register_interface(const std::string& interface_name,
                               std::function<std::vector<std::string>()> get_known_types,
                               std::function<Interface::pointer(const std::string&, bool)> instantiate)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            // Check if already registered
            for (const auto& info : m_interfaces) {
                if (info.interface_name == interface_name) {
                    return;  // Already registered
                }
            }
            m_interfaces.push_back({interface_name, get_known_types, instantiate});
        }

        /// Get all registered interface types
        std::vector<InterfaceInfo> all_interfaces() const
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_interfaces;
        }

        /// Singleton access
        static GlobalFactoryRegistry& instance()
        {
            static GlobalFactoryRegistry s_instance;
            return s_instance;
        }

       private:
        GlobalFactoryRegistry() {}
        std::vector<InterfaceInfo> m_interfaces;
        mutable std::mutex m_mutex;
    };

    /** A templated factory of objects of type Type that associates a
     * name to an object, returning a preexisting one if it exists. */
    template <class Type>
    class NamedFactory : public WireCell::INamedFactory {
       public:
        /// Remember the underlying type.
        typedef Type type;

        /// The exposed pointer type.
        typedef std::shared_ptr<Type> pointer_type;

        NamedFactory()
          : m_classname("")
        {
        }

        /// Return existing instance of given name or nullptr if not found.
        Interface::pointer find(const std::string& name)
        {
            auto it = m_objects.find(name);
            if (it == m_objects.end()) {
                return nullptr;
            }
            return it->second;
        }

        /// Return an instance associated with the given name.
        Interface::pointer create() { return create(""); }
        Interface::pointer create(const std::string& name)
        {
            auto it = m_objects.find(name);
            if (it == m_objects.end()) {
                pointer_type p(new Type);
                m_objects[name] = p;
                return p;
            }
            return it->second;
        }

        virtual void set_classname(const std::string& name) { m_classname = name; }
        virtual const std::string& classname() { return m_classname; }

       private:
        std::unordered_map<std::string, pointer_type> m_objects;
        std::string m_classname;
    };

    /** A registry of factories that produce instances which implement
     * a given interface. */
    template <class IType>
    class NamedFactoryRegistry {
        Log::logptr_t l;

       public:
        typedef IType interface_type;
        typedef std::shared_ptr<IType> interface_ptr;
        typedef WireCell::INamedFactory* factory_ptr;
        typedef std::unordered_map<std::string, factory_ptr> factory_lookup;
        typedef std::set<std::string> known_type_set;

        NamedFactoryRegistry()
          : l(Log::logger("factory"))
        {
            // Register this interface-specific registry with the global registry
            std::string interface_name = demangle(typeid(IType).name());

            // Capture this registry in lambdas for the global registry
            auto get_types = [this]() -> std::vector<std::string> {
                auto ktset = this->known_types();
                return std::vector<std::string>(ktset.begin(), ktset.end());
            };

            auto instantiate = [this](const std::string& classname, bool nullok) -> Interface::pointer {
                try {
                    interface_ptr iptr = this->instance(classname, "", true, nullok);
                    // Cast from interface_ptr (shared_ptr<IType>) to Interface::pointer (shared_ptr<Interface>)
                    return std::dynamic_pointer_cast<Interface>(iptr);
                }
                catch (...) {
                    if (nullok) return nullptr;
                    throw;
                }
            };

            GlobalFactoryRegistry::instance().register_interface(interface_name, get_types, instantiate);
        }
        size_t hello(const std::string& classname)
        {
            m_known_types.insert(classname);
            return m_known_types.size();
        }
        known_type_set known_types() const { return m_known_types; }

        /// Register an existing factory by the "class" name of the instance it can create.
        bool associate(const std::string& classname, factory_ptr factory)
        {
            m_lookup[classname] = factory;
            return true;
        }

        /// Look up an existing factory by the name of the "class" it can create.
        factory_ptr lookup_factory(const std::string& classname)
        {
            if (classname == "") {
                l->error("no class name given for type \"{}\"", demangle(typeid(IType).name()));
                return nullptr;
            }

            auto it = m_lookup.find(classname);
            if (it != m_lookup.end()) {
                return it->second;
            }

            // cache miss, try plugin

            WireCell::PluginManager& pm = WireCell::PluginManager::instance();

            std::string factory_maker = "make_" + classname + "_factory";
            auto plugin = pm.find(factory_maker.c_str());
            if (!plugin) {
                l->error("no plugin for \"{}\"", classname);
                return nullptr;
            }

            typedef void* (*maker_function)();
            maker_function mf;
            if (!plugin->symbol(factory_maker.c_str(), mf)) {
                l->error("no factory maker symbol for \"{}\"", classname);
                return nullptr;
            }

            void* fac_void_ptr = mf();

            if (!fac_void_ptr) {
                l->error("no factory for \"{}\"", classname);
                return nullptr;
            }

            factory_ptr fptr = reinterpret_cast<factory_ptr>(fac_void_ptr);
            m_lookup[classname] = fptr;
            return fptr;
        }

        /// Return instance of give type and optional instance name.
        /// If create is true, create the instance if it does not
        /// exist. If nullok is true return nullptr if it does not
        /// exist else throw by default.
        interface_ptr instance(const std::string& classname, const std::string& instname = "", bool create = true,
                               bool nullok = false)
        {
            factory_ptr fac = lookup_factory(classname);
            if (!fac) {
                if (nullok) {
                    return nullptr;
                }
                l->error("no factory for class \"{}\" (instance \"{}\")", classname, instname);
                std::cerr << "WireCell::NamedFactory: no factory for class \"" << classname << "\", (instance \""
                          << instname << "\"\n";

                THROW(FactoryException() << errmsg{"No factory for class " + classname});
            }
            WireCell::Interface::pointer iptr;
            std::string action = "";
            if (create) {
                iptr = fac->create(instname);
                action = "create";
            }
            else {
                iptr = fac->find(instname);
                action = "find";
            }
            if (!iptr) {
                if (nullok) {
                    return nullptr;
                }
                std::string msg = "NamedFactory: Failed to " + action + " instance \"" + instname;
                msg += "\" of class \"" + classname + "\"";
                l->error(msg);
                THROW(FactoryException() << errmsg{msg});
            }
            interface_ptr uptype = std::dynamic_pointer_cast<interface_type>(iptr);
            if (!uptype) {
                if (nullok) {
                    return nullptr;
                }
                std::string msg = "NamedFactory: Failed to cast instance: " + instname;
                msg += " of class " + classname;
                msg += " c++ type: " + type(iptr);
                msg += " to " + type(uptype);
                l->error(msg);
                THROW(FactoryException() << errmsg{msg});
            }
            return uptype;
        }

        /// Return a collection of class names known to this factory
        /// registry.  Note: linked/plugged shared libraries do not
        /// automatically register their factories.
        std::vector<std::string> known_classes()
        {
            std::vector<std::string> ret;
            for (auto it : m_lookup) {
                ret.push_back(it.first);
            }
            return ret;
        }

       private:
        factory_lookup m_lookup;
        known_type_set m_known_types;
    };

    /// Singleton interface
    namespace Factory {

        /// Associate a factory with the type it makes.
        template <class IType>
        bool associate(const std::string& classname, WireCell::INamedFactory* factory)
        {
            NamedFactoryRegistry<IType>& nfr = WireCell::Singleton<NamedFactoryRegistry<IType> >::Instance();
            bool ok = nfr.associate(classname, factory);
            if (ok) {
                return ok;
            }
            THROW(FactoryException() << errmsg{"Failed to associate class " + classname});
        }

        /// Lookup up a factory for the given type
        template <class IType>
        WireCell::INamedFactory* lookup_factory(const std::string& classname)
        {
            NamedFactoryRegistry<IType>& nfr = WireCell::Singleton<NamedFactoryRegistry<IType> >::Instance();
            WireCell::INamedFactory* ret = nfr.lookup_factory(classname);
            if (ret) {
                return ret;
            }
            THROW(FactoryException() << errmsg{"Failed to lookup factory for " + classname});
        }

        /// Lookup an interface by a type and optional name.  If
        /// create is true, create instance if missing.  if nullok is
        /// true return nullptr if can not create.  See also
        /// find<IType>().
        template <class IType>
        std::shared_ptr<IType> lookup(const std::string& classname, const std::string& instname = "",
                                      bool create = true, bool nullok = false)
        {
            NamedFactoryRegistry<IType>& nfr = WireCell::Singleton<NamedFactoryRegistry<IType> >::Instance();
            std::shared_ptr<IType> ret = nfr.instance(classname, instname, create, nullok);
            if (ret) {
                return ret;
            }
            if (nullok) {
                return nullptr;
            }
            THROW(FactoryException() << errmsg{"Failed to lookup instance for " + classname + " " + instname});
        }

        /// Return existing instance of given classname with instname
        /// Throws on failure.
        template <class IType>
        std::shared_ptr<IType> find(const std::string& classname, const std::string& instname = "")
        {
            std::shared_ptr<IType> ret = lookup<IType>(classname, instname, false, false);
            return ret;
        }
        /// As above but quietly returns nullptr on failure
        template <class IType>
        std::shared_ptr<IType> find_maybe(const std::string& classname, const std::string& instname = "")
        {
            std::shared_ptr<IType> ret = lookup<IType>(classname, instname, false, true);
            return ret;
        }

        /// Lookup an interface by a type:name pair.
        template <class IType>
        std::shared_ptr<IType> lookup_tn(const std::string& tn, bool create = true, bool nullok = false)
        {
            if (tn.empty()) {
                if (nullok) {
                    return nullptr;
                }
                THROW(FactoryException() << errmsg{"Empty type:name string"});
            }

            std::string t, n;
            std::tie(t, n) = String::parse_pair(tn);
            return lookup<IType>(t, n, create, nullok);
        }

        /// Like lookup_tn but with create false.
        template <class IType>
        std::shared_ptr<IType> find_tn(const std::string& tn)
        {
            std::shared_ptr<IType> ret = lookup_tn<IType>(tn, false, false);
            return ret;
        }
        /// Like find_tn but with nullok true
        template <class IType>
        std::shared_ptr<IType> find_maybe_tn(const std::string& tn)
        {
            std::shared_ptr<IType> ret = lookup_tn<IType>(tn, false, true);
            return ret;
        }

        /// Return a vector of all known classes of given interface.
        template <class IType>
        std::vector<std::string> known_classes()
        {
            NamedFactoryRegistry<IType>& nfr = WireCell::Singleton<NamedFactoryRegistry<IType> >::Instance();
            return nfr.known_classes();
        }
        template <class IType>
        std::vector<std::string> known_types()
        {
            NamedFactoryRegistry<IType>& nfr = WireCell::Singleton<NamedFactoryRegistry<IType> >::Instance();
            auto ktset = nfr.known_types();
            std::vector<std::string> ret(ktset.begin(), ktset.end());
            return ret;
        }

        /// Access the global factory registry for introspection
        inline GlobalFactoryRegistry& global_registry()
        {
            return GlobalFactoryRegistry::instance();
        }

        /// Get all registered interface types
        inline std::vector<GlobalFactoryRegistry::InterfaceInfo> all_interfaces()
        {
            return GlobalFactoryRegistry::instance().all_interfaces();
        }
    }  // namespace Factory

}  // namespace WireCell

template <class Concrete, class... Interface>
void* make_named_factory_factory(std::string name)
{
    static void* void_factory = nullptr;
    if (!void_factory) {
        void_factory = new WireCell::NamedFactory<Concrete>;
        WireCell::NamedFactory<Concrete>* factory = reinterpret_cast<WireCell::NamedFactory<Concrete>*>(void_factory);
        std::vector<bool> ret{WireCell::Factory::associate<Interface>(name, factory)...};
    }
    return void_factory;
}

template <class Concrete, class... Interface>
size_t namedfactory_hello(std::string name)
{
    std::vector<size_t> ret{WireCell::Singleton<WireCell::NamedFactoryRegistry<Interface> >::Instance().hello(name)...};
    return ret.size();
}

#define WIRECELL_FACTORY(NAME, CONCRETE, ...)                                                          \
    static size_t hello_##NAME##_me = namedfactory_hello<CONCRETE, __VA_ARGS__>(#NAME);                \
    extern "C" {                                                                                       \
    void* make_##NAME##_factory() { return make_named_factory_factory<CONCRETE, __VA_ARGS__>(#NAME); } \
    }

#endif
