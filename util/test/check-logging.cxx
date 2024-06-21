/** Excerpt of logging-related code from WireCell::Main

   This allows checking WCT logging as it would be used via wire-cell.

   Print help

   $ check-logging -h

   Default log stream, level and sink:

   $ check-logging

   Limit the default stream

   $ check-logging -L warn

   Setting any log sink (-l/--logsink) will cause the default log stream to be named "wct".

   Set a default log sink

   $ check-logging -l stderr

   Set a sinks with levels

   $ check-logging -l stderr:trace -l info.log:info

   Set message level

   $ check-logging -l stderr:trace -l info.log:info -l trace

   Set message level on per stream basis

   $ check-logging -l stderr:trace -l info.log:info -l trace -L wct:critical -L test:warn
   
*/

#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"

#include <boost/program_options.hpp>

#include <iostream>
#include <vector>
#include <string>

namespace po = boost::program_options;
using namespace WireCell;

static void add_logsink(const std::string& log, const std::string& level="")
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
static void set_loglevel(const std::string& log, const std::string& level="")
{
    // std::cerr << "Set log \""<<log<<"\" to level \"" << level << "\"\n";
    Log::set_level(level, log);
}


int main(int argc, char* argv[])
{
    // From Main::cmdline
    po::options_description desc("Test WCT logging\n\nUsage:\n\twire-cell [-h/--help|options] [message ...]\n\nOptions");
    desc.add_options()("help,h", "produce help message")

        ("logsink,l", po::value<std::vector<std::string> >(),
         "set log sink as <filename> or 'stdout' or 'stderr', "
         "a log level for the sink may be given by appending ':<level>'")

        ("loglevel,L", po::value<std::vector<std::string> >(),
         "set lowest log level for a log in form 'name:level' "
         "or just give 'level' value for all "
         "(level one of: critical,error,warn,info,debug,trace)")

        ;

    po::positional_options_description pos;
    pos.add("config", -1);

    po::variables_map opts;
    // po::store(po::parse_command_line(argc, argv, desc), opts);
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(pos).run(), opts);
    po::notify(opts);


    if (opts.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (opts.count("logsink")) {
        for (auto ls : opts["logsink"].as<std::vector<std::string> >()) {
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
        for (auto ll : opts["loglevel"].as<std::vector<std::string> >()) {
            auto lal = String::split(ll, ":");
            if (lal.size() == 2) {
                set_loglevel(lal[0], lal[1]);
            }
            else {
                set_loglevel("", lal[0]);
            }
        }
    }

    // From Main::initialize()
    set_loglevel("", "debug");

    {
        spdlog::level::level_enum active_level = (spdlog::level::level_enum)SPDLOG_ACTIVE_LEVEL;
        spdlog::info("compiled active level {} ({})",
                     spdlog::level::to_short_c_str(active_level),
                     SPDLOG_ACTIVE_LEVEL);
        spdlog::trace("default trace");
        spdlog::debug("default debug");
        SPDLOG_TRACE("default TRACE");
        SPDLOG_DEBUG("default DEBUG");
        spdlog::info("default info");
        spdlog::warn("default warn");
        spdlog::error("default error");
        spdlog::critical("default critical");

    }

    Log::logptr_t log(Log::logger("test", false)); // shared is false in Aux::Logger

    log->trace("test trace");
    log->debug("test debug");
    SPDLOG_LOGGER_TRACE(log,"test TRACE");
    SPDLOG_LOGGER_DEBUG(log,"test DEBUG");
    log->info("test info");
    log->warn("test warn");
    log->error("test error");
    log->critical("test critical");

    return 0;
}
