/**
   BlobDepoFill produces a new cluster with blob charge filled
   from depo charge.

   Except for replacing blob charge all input information is retained
   in the output.  Blobs in which no depo charge is found will get
   assigned zero value.  Depo charge not located inside any input
   blobs will be lost and in particular will not result in any new
   blobs.

   The the depos to input are best to be post-drift in order to more
   accurately align them to blob volumes and to filter them in the
   drift direction to be relevant to the anode planes under
   consideration.  OTOH, the depo set may span multiple AnodePlanes in
   the transverse direction as the input depos will be filtered so
   that only those in the active region of the anode plane are
   considered (ie, same as DepoTransform).  The faces referenced by
   the blobs are used to perform this filtering.

   See DepoBlobFill.org for details on how this component works.
 */

#ifndef WIRECELLIMG_BLOBDEPOFILL
#define WIRECELLIMG_BLOBDEPOFILL

#include "WireCellAux/Logger.h"

#include "WireCellIface/IBlobDepoFill.h"
#include "WireCellIface/IConfigurable.h"

#include <functional>

namespace WireCell::Img {

    class BlobDepoFill : public Aux::Logger,
                         public IBlobDepoFill,
                         public IConfigurable
    {
      public:

        BlobDepoFill();
        virtual ~BlobDepoFill();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // IBlobDepoFill.  Input tuple is (blob,depo).
        virtual bool operator()(const input_tuple_type& intup,
                                output_pointer& out);

      private:
        double m_toffset{0};
        double m_nsigma{3.0};
        double m_speed{0};      // must configure
        int m_pindex{2};        // primary plane index
        size_t m_count{0};
    };

}

#endif

