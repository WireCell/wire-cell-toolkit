/** Sample blobs to make point cloud tree and output as tensors.

    The heaving lifting is done via IBlobSamper instances.
*/
#ifndef WIRECELL_CLUS_BLOBSAMPLING
#define WIRECELL_CLUS_BLOBSAMPLING

#include "WireCellIface/IBlobSampling.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"


namespace WireCell::Clus {

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
        
        /** Configuration: "samplers"

            An object with attributes providing names of
            IBlobSamplers.  The attribute names will be used to name
            the point cloud produced by the samplers.
        */
        std::map<std::string, IBlobSampler::pointer> m_samplers;

        /** Config: "datapath"

            Set the datapath for the tensor representing the point
            cloud tree.  If a %d format code is found it wil be
            interpolated with the IBlobSet::ident() value.
         */
        std::string m_datapath = "pointtrees/%d";

        size_t m_count{0};

    };
}

#endif
