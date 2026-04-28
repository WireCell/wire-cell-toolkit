#include "WireCellImg/BlobLabeling.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(BlobLabeling, WireCell::Img::BlobLabeling,
                 WireCell::INamed,
                 WireCell::IFunctionNode<WireCell::IBlobSet, WireCell::IBlobSet>,
                 WireCell::IConfigurable)

using namespace WireCell;

Img::BlobLabeling::BlobLabeling()
    : Aux::Logger("BlobLabeling", "img")
{
}

Img::BlobLabeling::~BlobLabeling()
{
}

WireCell::Configuration Img::BlobLabeling::default_configuration() const
{
    Configuration cfg;
    // Add default configuration parameters here as needed
    return cfg;
}

void Img::BlobLabeling::configure(const WireCell::Configuration& cfg)
{
    // Configure the component with parameters from cfg
    // For now, this is empty as we have no configuration parameters
}

bool Img::BlobLabeling::operator()(const input_pointer& in, output_pointer& out)
{
    // Handle end-of-stream
    if (!in) {
        out = nullptr;
        log->debug("EOS");
        return true;
    }

    // For now, this is a simple pass-through
    // The input IBlobSet is directly assigned to output
    out = in;
    
    // log->debug("BlobLabeling: processed blob set with {} blobs", 
    //            in->blobs().size());

    return true;
}