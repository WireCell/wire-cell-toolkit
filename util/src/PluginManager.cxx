#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Exceptions.h"

#include <string>
#include <dlfcn.h>

using namespace WireCell;
using namespace std;

Plugin::Plugin(void* lib)
  : m_lib(lib)
{
}
Plugin::~Plugin() { dlclose(m_lib); }

void* Plugin::raw(const std::string& symbol_name)
{
    void* ret = dlsym(m_lib, symbol_name.c_str());
    return ret;
}

bool Plugin::contains(const std::string& symbol_name) { return nullptr != raw(symbol_name); }

PluginManager& WireCell::PluginManager::instance()
{
    static PluginManager inst;
    return inst;
}

WireCell::Plugin* WireCell::PluginManager::add(const std::string& plugin_name, const std::string& libname)
{
    Plugin* plugin = get(plugin_name);
    if (plugin) {
        l->debug("already have plugin {}", plugin_name);
        return plugin;
    }

    std::string exts[2] = {".so", ".dylib"};
    for (int ind = 0; ind < 2; ++ind) {
        std::string ext = exts[ind];
        string lname = "";
        if (libname == "") {
            lname = "lib";
            lname += plugin_name;
            lname += ext;
        }
        else {
            lname = libname;
        }
        // RTLD_NODELETE: keep the plugin mapped even across dlclose, so that
        // global factory makers (std::function targets living in this lib)
        // remain valid through static teardown at program exit. Without it,
        // the plugin is unmapped at exit before the GlobalFactoryRegistry
        // singleton is destroyed -> use-after-unload crash (macOS).
        void* lib = dlopen(lname.c_str(), RTLD_NOW | RTLD_NODELETE);
        if (!lib) {
            l->error("Failed to load {}: {}", lname, dlerror());
            continue;
        }

        m_plugins[plugin_name] = new Plugin(lib);
        l->debug("loaded plugin #{} \"{}\" from library \"{}\": {}", m_plugins.size(), plugin_name, lname,
                 (void*) m_plugins[plugin_name]);
        return m_plugins[plugin_name];
    }
    l->critical("failed to load plugin: \"{}\"", plugin_name);
    THROW(IOError() << errmsg{"failed to load plugin: " + plugin_name});
    return nullptr;
}

WireCell::Plugin* WireCell::PluginManager::get(const std::string& plugin_name)
{
    auto pit = m_plugins.find(plugin_name);
    if (pit == m_plugins.end()) {
        return nullptr;
    }
    return pit->second;
}

WireCell::Plugin* WireCell::PluginManager::find(const std::string& symbol_name)
{
    for (auto pit : m_plugins) {
        Plugin* maybe = pit.second;
        if (maybe->contains(symbol_name)) {
            return maybe;
        }
    }
    return nullptr;
}

WireCell::PluginManager::PluginManager()
  : l(Log::logger("sys"))
{
}
WireCell::PluginManager::~PluginManager()
{
    for (auto pit : m_plugins) {
        delete pit.second;
        pit.second = nullptr;
    }
}
