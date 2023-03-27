#ifndef WIRECELL_IMG_BLOBSETFANOUT
#define WIRECELL_IMG_BLOBSETFANOUT

#include "WireCellIface/IBlobSetFanout.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Img {

    // Fan out 1 blobset to N set at construction or configuration time.
    class BlobSetFanout : public Aux::Logger,
                          public IBlobSetFanout, public IConfigurable {
    public:
        BlobSetFanout(size_t multiplicity = 0);
        virtual ~BlobSetFanout();

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> output_types();

        // IFanout
        virtual bool operator()(const input_pointer& in, output_vector& outv);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    private:
        size_t m_multiplicity;
        size_t m_count{0};

            
    };
}

#endif
