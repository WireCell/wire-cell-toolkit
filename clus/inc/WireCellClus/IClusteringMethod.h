/**

   Interface providing a clustering method.

   Note, this interface is local to WCT's "clus/" sub-package and subject to
   change.

 */
#ifndef WIRECELLCLUS_ICLUSTERINGMETHOD
#define WIRECELLCLUS_ICLUSTERINGMETHOD


#include "WireCellUtil/IComponent.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

#include <set>

namespace WireCell::Clus {

    class IClusteringMethod : public IComponent<IClusteringMethod> {
    public:

        using cluster_set_t = std::set<const Facade::Cluster*>;
        
        virtual ~IClusteringMethod() {};

        /// Mutate a and/or b and maybe fill c.  Note, this API is subject to change
        virtual void clustering(Facade::Grouping& a, Facade::Grouping& b) const = 0;
    };
}


#endif
