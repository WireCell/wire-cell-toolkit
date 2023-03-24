// A wire cell version of the jsonnet cli.

#include "CLI11.hpp"

#include "WireCellUtil/Persist.h"
#include "WireCellUtil/String.h"

#include <iostream>
#include <vector>
#include <map>

using namespace WireCell;

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
    CLI::App app{"wcsonnet is a Wire-Cell Toolkit aware Jsonnet compiler"};

    std::string filename, output="/dev/stdout";
    std::vector<std::string> load_path, extvars, extcode, tlavars, tlacode;

    app.add_option("-o,--output", output,
                   "Output file");
    app.add_option("-P,--path", load_path,
                   "Search paths to consider in addition to those in WIRECELL_PATH")->type_size(1)->allow_extra_args(false);
    app.add_option("-V,--ext-str", extvars,
                   "Jsonnet external variable as <name>=<string>")->type_size(1)->allow_extra_args(false);
    app.add_option("-C,--ext-code", extcode,
                   "Jsonnet external code as <name>=<string>")->type_size(1)->allow_extra_args(false);
    app.add_option("-A,--tla-str", tlavars,
                   "Jsonnet level argument value <name>=<string>")->type_size(1)->allow_extra_args(false);
    app.add_option("-S,--tla-code", tlacode,
                   "Jsonnet level argument code <name>=<string>")->type_size(1)->allow_extra_args(false);
    app.add_option("file", filename, "Jsonnet file to compile");
                    
    // app.set_help_flag();
    // auto help = app.add_flag("-h,--help", "Print help message");

    CLI11_PARSE(app, argc, argv);

    // try {
    //     app.parse(argc, argv);
    //     if (*help) {
    //         throw CLI::CallForHelp();
    //     }
    // }
    // catch(const CLI::Error &e) {
    //     return app.exit(e);
    // }


    std::map<std::string, std::string> m_extvars, m_extcode, m_tlavars, m_tlacode;

    try {

        parse_param("--ext-str", extvars, m_extvars);
        parse_param("--ext-code", extcode, m_extcode);
        parse_param("--tla-str", tlavars, m_tlavars);
        parse_param("--tla-code", tlacode, m_tlacode);

        if (filename.empty()) {
            std::cerr << "Must give at least one Jsonnet file to compile" << std::endl;
            throw CLI::CallForHelp();
        }
    }
    catch(const CLI::Error &e) {
        return app.exit(e);
    }

    Persist::Parser parser(load_path, m_extvars, m_extcode,
                           m_tlavars, m_tlacode);

    auto jdat = parser.load(filename);
    std::ofstream out(output);
    out << jdat << std::endl;

    return 0;
}
