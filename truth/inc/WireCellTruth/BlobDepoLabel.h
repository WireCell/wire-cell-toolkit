/**
   BlobDepoLabel: BlobDepoFill with the per-track bookkeeping of
   DepoBlobAssoc, same interface (ICluster, IDepoSet) -> ICluster.

   With the extra smearing knobs off the output cluster is
   bit-identical to BlobDepoFill's, so it can replace BlobDepoFill in
   an existing truth tier (the "tru0"/"tru" catchers) under a byte
   gate.  The per-track table itself is emitted by BlobDepoLabelTensor;
   this node only logs its summary.

   wcp-porting-img pr_ml/docs/19 sec 3.5.
 */

#ifndef WIRECELLTRUTH_BLOBDEPOLABEL
#define WIRECELLTRUTH_BLOBDEPOLABEL

#include "WireCellTruth/DepoBlobAssoc.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IBlobDepoFill.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::Truth {

    class BlobDepoLabel : public Aux::Logger,
                          public IBlobDepoFill,
                          public IConfigurable
    {
      public:
        BlobDepoLabel();
        virtual ~BlobDepoLabel();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_tuple_type& intup,
                                output_pointer& out);

      private:
        DepoBlobAssocParams m_params;
        size_t m_count{0};
    };

    /// Shared by BlobDepoLabel and BlobDepoLabelTensor.
    Configuration assoc_default_configuration(const DepoBlobAssocParams& p);
    void assoc_configure(DepoBlobAssocParams& p, const Configuration& cfg);
}

#endif
