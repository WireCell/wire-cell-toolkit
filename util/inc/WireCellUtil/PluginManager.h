#ifndef WIRECELL_PLUGINMANAGER
#define WIRECELL_PLUGINMANAGER

#include "WireCellUtil/Logging.h"

#include <map>
#include <mutex>
#include <string>

namespace WireCell {

    class Plugin {
       public:
        Plugin(void* lib);
        ~Plugin();

        void* raw(const std::string& symbol_name);

        bool contains(const std::string& symbol_name);

        template <typename T>
        bool symbol(const std::string& symbol_name, T& ret)
        {
            void* thing = raw(symbol_name);
            if (!thing) {
                return false;
            }
            ret = reinterpret_cast<T>(thing);
            return true;
        }

       private:
        void* m_lib;
    };

    /** This is meant to be used from a WireCell::Singleton. */
    class PluginManager {
        Log::logptr_t l;
        PluginManager();
        ~PluginManager();

       public:
        static PluginManager& instance();

        /// Add a plugin.  If libname is not given, try to derive it.
        Plugin* add(const std::string& plugin_name, const std::string& libname = "");

        Plugin* get(const std::string& plugin_name);

        Plugin* find(const std::string& symbol_name);

       private:
        std::map<std::string, Plugin*> m_plugins;
        // Guards m_plugins: add/get/find may be called from several
        // threads (e.g. Factory lookups from configure() in more than one
        // WireCell::Main).  Recursive because add() calls get().
        mutable std::recursive_mutex m_mutex;
    };

}  // namespace WireCell

#endif
