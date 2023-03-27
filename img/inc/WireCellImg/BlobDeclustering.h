/** Forward blobs in cluster as a blob set.
*/
#ifndef WIRECELL_IMG_BLOBDECLUSTERING
#define WIRECELL_IMG_BLOBDECLUSTERING

#include "WireCellIface/IBlobDeclustering.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Img {

    class BlobDeclustering : public Aux::Logger, public IBlobDeclustering
    {
      public:
        BlobDeclustering();
        virtual ~BlobDeclustering();

        virtual bool operator()(const input_pointer& blobset, output_pointer& tensorset);

      private:
        
        size_t m_count{0};
    };
}

#endif
