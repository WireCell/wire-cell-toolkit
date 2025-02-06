#ifndef WIRECELL_LOGGING
#define WIRECELL_LOGGING

// See this header for info about build configure and compile time control over
// SPDLOG levels.
#include "WireCellUtil/Spdlog.h"

#include <memory>
#include <string>

namespace WireCell {

    namespace Log {

        typedef std::shared_ptr<spdlog::logger> logptr_t;
        typedef std::shared_ptr<spdlog::sinks::sink> sinkptr_t;

        // WCT maintains a shared collection of sinks associated with
        // its default loggers created through this API.  No sinks are
        // added by default.  The WCT application should add some if
        // output is wanted.  All sinks added will be applied to any
        // subsequently made loggers by logger() below.

        // Add a log file sink with optional level.
        void add_file(std::string filename, std::string level = "");

        // Add a standard out console sink with optional level.
        void add_stdout(bool color = true, std::string level = "");

        // Add a standard err console sink with optional level.
        void add_stderr(bool color = true, std::string level = "");

        // Return a logger by name, making it if it does not yet
        // exist.  If shared_sinks is true, the logger will be
        // attached to a shared set of sinks created by prior calls to
        // the above add_*() functions.  If false then copies of the
        // previously added sinks will be attached (which will be
        // needed if custom patterns will be set on the logger).  WCT
        // components are encouraged to may make unique loggers with
        // some short name related to the component type/name and hold
        // on to them for use in their code.
        logptr_t logger(std::string name, bool share_sinks=true);

        // Set log level.  If which is empty the set level of logs.
        // Otherwise, set the given logger.
        void set_level(std::string level, std::string which = "");

        // W.r.t. levels, loggers can have their level explicitly set with
        // set_level() or implicitly set (to default level).  Calling this
        // function will explicitly set the level for an otherwise implicitly
        // set logger to have the level of the best matching explicitly set
        // logger.  Best matching is defined as the longest name of an
        // explicitly set logger that is the prefix of the name of the
        // implicitly set logger.  Note, calling this does not cause an
        // implicitly set logger to become explicitly set.
        void fill_levels();

        // Set logging pattern the default or given logger's sinks.
        void set_pattern(std::string pattern, std::string which = "");

        // Set up logging.  This is intended for a oneliner in main().  Output
        // can be "stderr", "stdout" or a file name.
        void default_logging(const std::string& output = "stderr",
                             std::string level = "", bool with_env = true);

    }  // namespace Log

}  // namespace WireCell

#endif

