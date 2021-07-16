/** Single-point main entry to Wire Cell Toolkit.
 */

#include "WireCellApps/Main.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Point.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IApplication.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <vector>
#include <iostream>

#include "config.h"

#if HAVE_FFTWTHREADS_LIB
#include <fftw3.h>
#endif

#if HAVE_TBB_LIB
#include <tbb/global_control.h>
#endif

using namespace WireCell;
using namespace std;
namespace po = boost::program_options;
using namespace boost::algorithm;
using namespace boost::property_tree;

Main::Main()
  : l(Log::logger("main"))
{
#if HAVE_FFTWTHREADS_LIB
    fftwf_init_threads();
    fftwf_make_planner_thread_safe();
#endif
}

Main::~Main() { finalize(); }

int Main::cmdline(int argc, char* argv[])
{
    // clang-format off
    po::options_description desc("Options");
    desc.add_options()("help,h", "wire-cell [options] [arguments]")

        ("logsink,l", po::value<vector<string> >(),
         "set log sink as <filename> or 'stdout' or 'stderr', "
         "a log level for the sink may be given by appending ':<level>'")

        ("loglevel,L", po::value<vector<string> >(),
         "set lowest log level for a log in form 'name:level' "
         "or just give 'level' value for all "
         "(level one of: critical,error,warn,info,debug,trace)")

        ("app,a", po::value<vector<string> >(),
         "application component to invoke")

        ("config,c", po::value<vector<string> >(),
         "provide a configuration file")
        
        ("plugin,p", po::value<vector<string> >(),
         "specify a plugin as name[:lib]")
        
        ("ext-str,V", po::value<vector<string> >(),
         "specify a Jsonnet external variable=<string>")
        
        ("ext-code,C", po::value<vector<string> >(),
         "specify a Jsonnet external variable=<code>")
        
        ("tla-str,A", po::value<vector<string> >(),
         "specify a Jsonnet top level arguments variable=<string>")
        
        ("tla-code", po::value<vector<string> >(),
         "specify a Jsonnet top level arguments variable=<code>")
        
        ("path,P", po::value<vector<string> >(),
         "add to JSON/Jsonnet search path")
#if HAVE_TBB_LIB
        ("threads,t", po::value<int>(),
         "limit number of threads used")
#endif
        ;
    // clang-format on

    po::variables_map opts;
    po::store(po::parse_command_line(argc, argv, desc), opts);
    po::notify(opts);

    if (opts.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (opts.count("config")) {
        for (auto fname : opts["config"].as<vector<string> >()) {
            add_config(fname);
        }
    }

    if (opts.count("path")) {
        for (auto path : opts["path"].as<vector<string> >()) {
            add_path(path);
        }
    }

    // Get any external variables
    if (opts.count("ext-str")) {
        for (auto vev : opts["ext-str"].as<vector<string> >()) {
            auto vv = String::split(vev, "=");
            add_var(vv[0], vv[1]);
        }
    }
    // And any external code
    if (opts.count("ext-code")) {
        for (auto vev : opts["ext-code"].as<vector<string> >()) {
            auto vv = String::split(vev, "=");
            add_code(vv[0], vv[1]);
        }
    }
    // Add any top-level argument string values
    if (opts.count("tla-str")) {
        for (auto vev : opts["tla-str"].as<vector<string> >()) {
            auto vv = String::split(vev, "=");
            tla_var(vv[0], vv[1]);
        }
    }
    // Add any top-level argument code
    if (opts.count("tla-code")) {
        for (auto vev : opts["tla-code"].as<vector<string> >()) {
            auto vv = String::split(vev, "=");
            tla_code(vv[0], vv[1]);
        }
    }

    // fixme: these aren't yet supported.
    // if (opts.count("jsonpath")) {
    //     jsonpath_vars = opts["jsonpath"].as< vector<string> >();
    // }

    if (opts.count("plugin")) {
        for (auto plugin : opts["plugin"].as<vector<string> >()) {
            add_plugin(plugin);
        }
    }
    if (opts.count("app")) {
        for (auto app : opts["app"].as<vector<string> >()) {
            add_app(app);
        }
    }
    if (opts.count("logsink")) {
        for (auto ls : opts["logsink"].as<vector<string> >()) {
            auto ll = String::split(ls, ":");
            if (ll.size() == 1) {
                add_logsink(ll[0]);
            }
            if (ll.size() == 2) {
                add_logsink(ll[0], ll[1]);
            }
        }
    }

    if (opts.count("loglevel")) {
        for (auto ll : opts["loglevel"].as<vector<string> >()) {
            auto lal = String::split(ll, ":");
            if (lal.size() == 2) {
                set_loglevel(lal[0], lal[1]);
            }
            else {
                set_loglevel("", lal[0]);
            }
        }
    }

#ifdef HAVE_TBB_LIB
    if (opts.count("threads")) {
        m_threads = opts["threads"].as<int>();
    }
#endif

    // Maybe make this cmdline configurable.  For now, set all
    // backends the same.
    Log::set_pattern("[%H:%M:%S.%03e] %L [%^%=8n%$] %v");

    return 0;
}

void Main::add_plugin(const std::string& libname) { m_plugins.push_back(libname); }

void Main::add_app(const std::string& tn) { m_apps.push_back(tn); }

void Main::add_logsink(const std::string& log, const std::string& level)
{
    if (log == "stdout") {
        Log::add_stdout(true, level);
        return;
    }
    if (log == "stderr") {
        Log::add_stderr(true, level);
        return;
    }
    Log::add_file(log, level);
}
void Main::set_loglevel(const std::string& log, const std::string& level) { Log::set_level(level, log); }
void Main::add_config(const std::string& filename) { m_cfgfiles.push_back(filename); }

void Main::add_var(const std::string& name, const std::string& value) { m_extvars[name] = value; }

void Main::add_code(const std::string& name, const std::string& value) { m_extcode[name] = value; }

void Main::tla_var(const std::string& name, const std::string& value) { m_tlavars[name] = value; }

void Main::tla_code(const std::string& name, const std::string& value) { m_tlacode[name] = value; }

void Main::add_path(const std::string& dirname) { m_load_path.push_back(dirname); }

void Main::initialize()
{
    for (auto filename : m_cfgfiles) {
        l->info("loading config file {} ...", filename);
        Persist::Parser p(m_load_path, m_extvars, m_extcode, m_tlavars, m_tlacode);
        Json::Value one = p.load(filename);  // throws
        m_cfgmgr.extend(one);
        l->info("...done");
    }

    // Find if we have our special configuration entry
    int ind = m_cfgmgr.index("wire-cell");
    Configuration main_cfg = m_cfgmgr.pop(ind);
    if (!main_cfg.isNull()) {
        for (auto plugin : get<vector<string> >(main_cfg, "data.plugins")) {
            l->info("config requests plugin: \"{}\"", plugin);
            m_plugins.push_back(plugin);
        }
        for (auto app : get<vector<string> >(main_cfg, "data.apps")) {
            l->info("config requests app: \"{}\"", app);
            m_apps.push_back(app);
        }
    }

    // Load any plugin shared libraries requested by user.
    PluginManager& pm = PluginManager::instance();
    for (auto plugin : m_plugins) {
        string pname, lname;
        std::tie(pname, lname) = String::parse_pair(plugin);
        l->info("adding plugin: \"{}\"", plugin);
        if (lname.size()) {
            l->info("\t from library \"{}\"", lname);
        }
        pm.add(pname, lname);
    }

    // Apply any user configuration.  This is a two step.  First, just
    // assure all the components referenced in the configuration
    // sequence can be instantiated.  Then, find them again and
    // actually configure them.  This way, any problems fails fast.

    for (auto c : m_cfgmgr.all()) {
        if (c.isNull()) {
            continue;  // allow and ignore any totally empty configurations
        }
        if (c["type"].isNull()) {
            l->critical("all configuration must have a type attribute, got: {}", c);
            THROW(ValueError() << errmsg{"got configuration sequence element lacking a type"});
        }
        string type = get<string>(c, "type");
        string name = get<string>(c, "name");
        l->info("constructing component: \"{}\":\"{}\"", type, name);
        auto iface = Factory::lookup<Interface>(type, name);  // throws
    }
    for (auto c : m_apps) {
        l->info("constructing app: \"{}\"", c);
        Factory::lookup_tn<IApplication>(c);
    }
    for (auto c : m_cfgmgr.all()) {
        if (c.isNull()) {
            continue;  // allow and ignore any totally empty configurations
        }
        string type = get<string>(c, "type");
        string name = get<string>(c, "name");
        l->info("configuring component: \"{}\":\"{}\"", type, name);
        auto cfgobj = Factory::find_maybe<IConfigurable>(type, name);  // doesn't throw.
        if (!cfgobj) {
            continue;
        }

        // Get component's hard-coded default config, update it with
        // anything the user may have provided and apply it.
        Configuration cfg = cfgobj->default_configuration();
        cfg = update(cfg, c["data"]);
        cfgobj->configure(cfg);  // throws
    }
}

void Main::operator()()
{
    // Find all IApplications to execute
    vector<IApplication::pointer> app_objs;
    for (auto component : m_apps) {
        string type, name;
        std::tie(type, name) = String::parse_pair(component);
        auto a = Factory::find<IApplication>(type, name);  // throws
        app_objs.push_back(a);
    }
    l->debug("executing {} apps, thread limit {}:",
             m_apps.size(), m_threads);

#if HAVE_TBB_LIB
    std::unique_ptr<tbb::global_control> gc;
    if (m_threads) {
        gc = std::make_unique<tbb::global_control>(
            tbb::global_control::max_allowed_parallelism,
            m_threads);
    }
    l->debug("executing {} apps, thread limit {}:",
             m_apps.size(), m_threads);
#else
    l->debug("executing {} apps:", m_apps.size());
#endif

    for (size_t ind = 0; ind < m_apps.size(); ++ind) {
        auto aobj = app_objs[ind];
        l->debug("executing app: \"{}\"", m_apps[ind]);
        aobj->execute();  // throws
    }

}

void Main::finalize()
{
    for (auto c : m_cfgmgr.all()) {
        if (c.isNull()) {
            continue;  // allow and ignore any totally empty configurations
        }
        string type = get<string>(c, "type");
        string name = get<string>(c, "name");
        auto doomed = Factory::find_maybe<ITerminal>(type, name);  // doesn't throw.
        if (!doomed) {
            continue;
        }
        l->info("finalizing component: \"{}\":\"{}\"", type, name);

        doomed->finalize();
    }
}
