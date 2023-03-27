// A wire cell CLI to operate on description of wires

#include "CLI11.hpp"

#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Exceptions.h"

#include <iostream>
#include <vector>
#include <map>

using namespace WireCell;
using namespace WireCell::WireSchema;

static void
parse_param(std::string name,
            const std::vector<std::string>& args,
            std::map<std::string,std::string>& store)
{
    for (auto one : args) {
        auto two = String::split(one, "=");
        if (two.size() != 2) {
            std::cerr
                << name
                << ": parameters are set as <name>=<value>, got "
                << one << std::endl;
            throw CLI::CallForHelp();
        }
        store[two[0]] = two[1];
    }
}

int main(int argc, char** argv)
{
    CLI::App app{"wcwires converts and validates Wire-Cell Toolkit wire descriptions"};

    std::string filename="", output="";
    std::vector<std::string> load_path;
    int correction = static_cast<int>(Correction::pitch);
    bool do_validate = false;
    bool do_fail_fast = false;
    double repsilon = 1e-6;

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

    if (output.empty()) {
        if (do_validate) {
            Store raw = load(filename.c_str(), Correction::load);
            try {
                validate(raw, repsilon, do_fail_fast);
            }
            catch (...) {
                std::cerr << "input invalid\n";
                return 1;
            }
            std::cerr << "input valid\n";
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
        std::cerr << "Loading " << filename << " with " << corstr << " corrections\n";

        Store store = load(filename.c_str(), cor);
        if (do_validate) {
            try {
                validate(store, repsilon, do_fail_fast);
            }
            catch (...) {
                std::cerr << "output invalid\n";
                return 1;
            }
            std::cerr << "output valid\n";
        }
        dump(output.c_str(), store);
    }
    return 0;
}
