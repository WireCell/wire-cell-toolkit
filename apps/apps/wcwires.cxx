// A wire cell CLI to operate on description of wires

#include "CLI11.hpp"

#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Logging.h"

#include <vector>
#include <string>
#include <cstdlib>              // setenv

using namespace WireCell;
using namespace WireCell::WireSchema;

using spdlog::debug;
using spdlog::info;
using spdlog::error;

int main(int argc, char** argv)
{
    CLI::App app{"wcwires converts and validates Wire-Cell Toolkit wire descriptions"};

    std::string filename="", output="";
    std::vector<std::string> load_path;
    int correction = static_cast<int>(Correction::pitch);
    bool do_validate = false;
    bool do_fail_fast = false;
    double repsilon = 1e-6;

    std::string logsink = "stderr";
    std::string loglevel = "info";

    app.add_option("-l,--logsink", logsink,
                   "Set log sink as <filename> or 'stdout' or 'stderr' (default)"
        )->type_size(1)->allow_extra_args(false);

    app.add_option("-L,--loglevel", loglevel,
                   "Set log level as 'debug' or 'info' (default)"
        )->type_size(1)->allow_extra_args(false);

    app.add_option("-P,--path", load_path,
                   "Search paths to consider in addition to those in WIRECELL_PATH"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-o,--output", output,
                   "Write out a wires file (def=none)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-c,--correction", correction,
                   "Correction level: 1=load,2=order,3=direction,4=pitch (def=4)"
        )->type_size(1)->allow_extra_args(false);
    app.add_flag("-v,--validate", do_validate,
                 "Perform input validation (def=false)"
        )->type_size(1)->allow_extra_args(false);
    app.add_flag("-f,--fail-fast", do_fail_fast,
                 "Fail on first validation error (def=false)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-e,--epsilon", repsilon,
                 "Unitless relative error determining imprecision during validation (def=1e-6)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("file", filename, "wires file");
                    
    // app.set_help_flag();
    // auto help = app.add_flag("-h,--help", "Print help message");

    CLI11_PARSE(app, argc, argv);

    Log::default_logging(logsink, loglevel, true);
    Log::set_level(loglevel);
    debug("logging to {} at level {}", logsink, loglevel);

    for (const auto& path : load_path) {
        std::string wpath = getenv("WIRECELL_PATH");
        debug("adding to WIRECELL_PATH: {}", path);
        if (wpath.empty()) {
            wpath = path;
        }
        else {
            wpath = path + ":" + wpath;
        }
        setenv ("WIRECELL_PATH", path.c_str(), 1);
    }


    if (output.empty()) {
        if (do_validate) {
            Store raw = load(filename.c_str(), Correction::load);
            try {
                validate(raw, repsilon, do_fail_fast);
            }
            catch (...) {
                error("wires file {} is invalid", filename);
                return 1;
            }
            info("valid");
        }
    }
    else {
        std::string corstr="";
        auto cor = static_cast<Correction>(correction);
        switch(cor) {
            case Correction::empty:
                corstr="load";
                cor = Correction::load;
                break;
            case Correction::load:
                corstr="load";
                break;
            case Correction::order:
                corstr="order";
                break;
            case Correction::direction:
                corstr="direction";
                break;
            case Correction::pitch:
                corstr="pitch";
                break;
            default:
                corstr="pitch";
                cor = Correction::pitch;
                break;
        };
        debug("Loading {} with {} corrections", filename, corstr);

        Store store = load(filename.c_str(), cor);
        if (do_validate) {
            try {
                validate(store, repsilon, do_fail_fast);
            }
            catch (...) {
                error("wires file {} is invalid", filename);
                return 1;
            }
            info("valid");
        }
        dump(output.c_str(), store);
    }
    return 0;
}
