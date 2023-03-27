/** Sample blobs to make point cloud datasets as tensor sets.

 The heaving lifting is done via an IBlobSamper.
*/
#ifndef WIRECELL_IMG_BLOBSAMPLING
#define WIRECELL_IMG_BLOBSAMPLING

#include "WireCellIface/IBlobSampling.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"


namespace WireCell::Img {

    class BlobSampling : public Aux::Logger, public IBlobSampling, public IConfigurable
    {
      public:
        BlobSampling();
        virtual ~BlobSampling();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& blobset, output_pointer& tensorset);

      private:
        
        /** Configuration: "sampler"

            Set the blob samper.
        */
        IBlobSampler::pointer m_sampler;

        /** Config: "datapath"

            Set the datapath for the tensor representing the frame.
            The string may provide a %d format code which will be
            interpolated with the frame's ident number.
         */
        std::string m_datapath = "frames/%d";

        size_t m_count{0};
    };
}

#endif
