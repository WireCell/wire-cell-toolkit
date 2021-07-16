/** This is the main entry point to the WCT.  It provides a single
 * command line method or optionally a family of fine-grained methods
 * for setup and running.
 *
 * One use of this is in the `wire-cell` command line program.  It may
 * also be use to embed WCT into a larger application or framework. */

#ifndef WIRECELL_MAIN
#define WIRECELL_MAIN

#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Persist.h"

#include "WireCellUtil/Logging.h"

#include <string>
#include <vector>

namespace WireCell {
    class Main {
       public:
        Main();
        ~Main();

        /// Single-point entry to Wire Cell.
        ///
        /// Pass in literal command line arguments and return a return
        /// code.  See `wire-cell --help` for args.
        ///
        /// Or, one can use subsequent methods for more fine-grained
        /// setup and execution.
        int cmdline(int argc, char* argv[]);

        /// Individual setup methods called by cmdline() or called
        /// explicitly by external application/frameork:

        /// Add an IApplication component to execute as a `type:name`
        /// string.
        void add_app(const std::string& tn);

        /// Append a top-level JSON/Jsonnet configuration file to the
        /// configuration sequence.
        void add_config(const std::string& filename);

        /// Bind an external scalar value to a variable so that it may
        /// be referenced in the configuration files (via
        /// std.extVar()).
        void add_var(const std::string& name, const std::string& value);

        /// Bind external configuration code (in Jsonnet language) to
        /// a variable so that its may be referenced the configuration
        /// files (via std.extVar())
        void add_code(const std::string& name, const std::string& value);

        /// Bind a "top level argument" to a string.  It only makes
        /// sense to use this if the main Jsonnet file results in a
        /// function with this argument.
        void tla_var(const std::string& name, const std::string& value);

        /// Bind a "top level argument" to some Jsonnet code.  It only
        /// makes sense to use this if the main Jsonnet file results
        /// in a function with this argument.
        void tla_code(const std::string& name, const std::string& value);

        // fixme: "jsonnet" CLI lets TLA code and strings also given in files....

        /// Add an element to the configuration path in which
        /// configuration files may be found.
        void add_path(const std::string& dirname);

        /// Add a plugin library in which components may be found.
        /// The libname may lack the initial "lib" prefix and file
        /// extension.
        void add_plugin(const std::string& libname);

        /// Add a log sink, reserved names 'stdout' and 'stderr' or a filename.
        void add_logsink(const std::string& log, const std::string& level = "");

        /// Set a minimum level to emit a message for a given
        /// log. (levels: critical, error, warn, info, debug, trace).
        void set_loglevel(const std::string& log, const std::string& level = "");

        /// Call once after all setup has been done and before
        /// running.
        void initialize();

        /// Run any and all application components once.
        void operator()();

        /// Call to finalize any terminal components.
        void finalize();



       private:
        ConfigManager m_cfgmgr;
        std::vector<std::string> m_plugins, m_apps, m_cfgfiles, m_load_path;
        Persist::externalvars_t m_extvars, m_extcode, m_tlavars, m_tlacode;
        Log::logptr_t l;

        // Limit number of threads.  0 means set no limit.  This is
        // only relevant if we are built with TBB support
        int m_threads{0};
    };

}  // namespace WireCell
#endif
