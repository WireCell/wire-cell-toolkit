/**
   These hold steiner-related free functions.

 */

#ifndef WIRECELLCLUS_STEINERFUNCTIONS
#define WIRECELLCLUS_STEINERFUNCTIONS

#include "SteinerGrapher.h"

namespace WireCell::Clus::Steiner {

    // Xin, these certainly require the "..."'s to be expanded as we figure out
    // WCT data types that are equivalent to what WCP uses.

    // Xin, in keeping with "grapher" == "pr3dcluster" I name these accordingly.
    // If they do NOT need the extra part of "grapher" it would be better to
    // pass the underlying Cluster and then rename the functions accordingly.

    void improve_grapher_2(Grapher& grapher/*, ...*/);
    void improve_grapher_1(Grapher& grapher/*,...*/);
    void improve_grapher(Grapher& grapher/*,...*/);
    void improve_grapher(Grapher& grapher, Grapher& other_grapher/*,...*/);


}


#endif 
