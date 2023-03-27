/**
   Some helper functions operating on IDepo
*/
#ifndef WIRECELLAUX_DEPOTOOLS
#define WIRECELLAUX_DEPOTOOLS

#include "WireCellIface/IDepo.h"
#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/Array.h"
#include "WireCellUtil/Binning.h"

namespace WireCell::Aux {

    /** Fill data and info arrays from IDepoSet.
             
        The "data" array will be shaped (ndepos,7) and each
        7-tuple holds: (time, charge, x, y, z, long, tran).

        The "info" array is integer and shaped (nedpos, 4).  Each
        4-tuple holds: (id, pdg, gen, child).

        If "priors" is true, then run the depos through flatten() in
        order to save them and give non-zero gen and child.  O.w.,
        only those depos in pass directly in the vector are saved.
    */
    void fill(Array::array_xxf& data, Array::array_xxi& info, 
              const IDepo::vector& depos, bool priors = true);
         

    /** Sort depos into the time bins that their Gaussian extent cover.
        
        Any excursion into a bin will cause the depo to be added.
     */ 
    std::vector<IDepo::vector>
    depos_by_slice(const IDepo::vector& depos, const Binning& tbins,
                   double time_offset=0, double nsigma=3.0);

    /** Return subset of depos which are within the sensitive area of
     * a single anode face or both faces of an anode plane.
     */
    IDepo::vector sensitive(const IDepo::vector& depos, IAnodeFace::pointer face);
    IDepo::vector sensitive(const IDepo::vector& depos, IAnodePlane::pointer anode);

}

#endif
