#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Testing.h"

using namespace WireCell;

int main(int argc, char* argv[])
{
    std::string fname = "microboone-celltree-wires-v2.1.json.bz2";
    int clevel = (int)WireSchema::Correction::pitch;
    if (argc > 1) {
        fname = argv[1];
    }
    if (argc > 2) {
        clevel = atoi(argv[2]);
    }

    if (clevel <= (int)WireSchema::Correction::empty) {
        clevel = (int)WireSchema::Correction::load;
    }
    if (clevel >= (int)WireSchema::Correction::ncorrections) {
        clevel = (int)WireSchema::Correction::pitch;
    }
    WireSchema::Correction cor = (WireSchema::Correction)clevel;

    std::cerr << "loading "<<fname<<" with correction " << (int)cor << ".\n";
    auto store = WireSchema::load(fname.c_str(), cor);
    try {
        WireSchema::validate(store);
    }
    catch (ValueError& err) {
        std::cerr << "invalid: "<<fname<<" with correction " << (int)cor << ".\n";
        //std::cerr << err.what() << std::endl;
        return -1;
    }
    std::cerr << "valid: "<<fname<<" with correction " << (int)cor << ".\n";
    return 0;
}
