/**

   Interface providing a clustering method.

   Note, this interface is local to WCT's "clus/" sub-package and subject to
   change.

 */
#ifndef WIRECELLCLUS_ICLUSTERINGMETHOD
#define WIRECELLCLUS_ICLUSTERINGMETHOD


#include "WireCellUtil/IComponent.h"
#include "WireCellClus/Facade_Ensemble.h"

#include <set>

namespace WireCell::Clus {

    class IClusteringMethod : public IComponent<IClusteringMethod> {
    public:

        virtual ~IClusteringMethod() {};

        /// Mutate the ensemble of cluster groupings.
        virtual void clustering(Facade::Ensemble& ensemble) const = 0;
    };
}


#endif
