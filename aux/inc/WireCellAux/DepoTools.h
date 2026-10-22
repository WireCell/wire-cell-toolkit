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

        Row order with priors: an active depo (gen 0) is followed
        immediately by its prior (gen 1), that one's prior (gen 2),
        etc.  The "child" of a gen>0 row is the row index of the depo
        whose prior() it is.  (Files written before 2026-10-21 have
        child == own row index for every gen>0 row; readers detect
        that pattern and fall back to "child = previous row".)
    */
    void fill(Array::array_xxf& data, Array::array_xxi& info,
              const IDepo::vector& depos, bool priors = true);

    /** Fill the optional "truth" sidecar array, row-aligned with
        the arrays of fill() (same flattening).

        The array is double and shaped (ndepos, 3).  Each 3-tuple
        holds: (energy, orig_id, reserved).  IDepo carries no
        pre-merge id, so orig_id and reserved are written as 0 and
        are reserved for a Stage-A side table.
    */
    void fill_truth(Eigen::ArrayXXd& truth,
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
