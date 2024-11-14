/*
  This provides a 2D version of an interval tree.

  For concepts see:

  https://en.wikipedia.org/wiki/Interval_tree

  It is based on Boost's ICL which provides 1D interval trees:

  https://www.boost.org/doc/libs/1_80_0/libs/icl/doc/html/index.html

  In short, it lets you do this:

                                              
  ┌──────┐                  ┌──────┐
  │ a    │                  └──────┘a
  │  ┌───┼───┐              ┌─┐┌───┐┌──┐
  │  │   │   │   disjoin    │a││a,b││b │
  └──┼───┘   │  ─────────>  └─┘└───┘└──┘
     │ b     │                 ┌───────┐
     └───────┘                 └───────┘b


  Thanks to herbstluftwm(1) for the diagram.

 */

#ifndef WIRECELL_UTIL_RECTANGLES
#define WIRECELL_UTIL_RECTANGLES

#include "WireCellUtil/Binning.h"

#if defined __GNUC__
 #pragma GCC diagnostic push
 #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#if defined __clang__
 #pragma clang diagnostic push
 #pragma clang diagnostic ignored "-Wunused-local-typedefs"
#endif

#include "boost/icl/interval_map.hpp"

#if defined __GNUC__
 #pragma GCC diagnostic pop
#endif

#if defined __clang__
 #pragma clang diagnostic pop
#endif

#include <vector>
#include <set>

namespace WireCell {

    /*
      Track rectangular "regions" (2D intervals).

      Key gives the type for rectangular coordinates.

      Value gives the type of data associated to regions.

      Set is used to collect Value objects and needs operator+= defined.

      See boost::icl docs for details.

     */
    // fixme: could be extended to N dimensions with template<int N>
    // and lots of variadic templating.
    template<typename Value,
             typename XKey,
             typename YKey = XKey,
             typename Set = std::set<Value>>
    class Rectangles {
      public:
        using xkey_t = XKey;
        using ykey_t = YKey;
        using value_t = Value;
        using set_t = Set;
        template<typename Key>
        using interval_t = typename boost::icl::interval<Key>::interval_type;
        using xinterval_t = interval_t<xkey_t>;
        using yinterval_t = interval_t<ykey_t>;

        // Intervals along the vertical Y-axis
        using ymap_t = boost::icl::interval_map<ykey_t, set_t>;
        // Intervals along the horizontal X-axis
        using xmap_t = boost::icl::interval_map<xkey_t, ymap_t>;
        // A rectangle is a pair of intervals: (xi,yi).
        using rectangle_t = std::pair<xinterval_t, yinterval_t>;
        // An element in the map is a pair of rectangle and value
        using element_t = std::pair<rectangle_t, value_t>;
        // A region is a contiguous rectangle and the set of values
        using region_t = std::pair<rectangle_t, set_t>;

        // Add a half-open rectangle.  Eg, width from [x1,x2), etc height.
        void add(const xkey_t& x1, const xkey_t& x2,
                 const ykey_t& y1, const ykey_t& y2,
                 const value_t& val)
        {
            add(xinterval_t::right_open(x1,x2),
                yinterval_t::right_open(y1,y2), val);
        }

        // Add rectangle as pair of intervals.
        void add(const xinterval_t& xi, const yinterval_t& yi, const set_t& vset) {
            ymap_t ymap;
            ymap += std::make_pair(yi, vset);
            m_xmap += std::make_pair(xi, ymap);
        }
        void add(const xinterval_t& xi, const yinterval_t& yi, const value_t& val) {
            add(xi, yi, set_t{val});
        }

        // Add rectangle as rectangle.
        void add(const rectangle_t& r, const set_t& vset) {
            add(r.first, r.second, vset);
        }
        void add(const rectangle_t& r, const value_t& val) {
            add(r.first, r.second, val);
        }

        // icl-like interface
        void operator+=(const element_t& ele) {
            add(ele.first, ele.second);
        }
        
        // Access the top interval map along X-axis.
        const xmap_t& xmap() const { return m_xmap; }

        // Return all contiguous regions of the overlapping
        // rectangles.
        std::vector<region_t> regions() const 
        {
            std::vector<region_t> ret;
            for (const auto& [xi, ymap] : m_xmap) {
                for (const auto& [yi, cs] : ymap) {
                    ret.push_back(region_t(rectangle_t(xi, yi), cs));
                }
            }
            return ret;
        }

        // Return regions that are the intersectection stored
        // rectangles and the given rectangle.
        std::vector<region_t> intersection(const rectangle_t& rec) const
        {
            return intersection(rec.first, rec.second);
        }
        // As above but for x/y intervals.
        std::vector<region_t> intersection(const xinterval_t& xi, const yinterval_t& yi) const
        {
            std::vector<region_t> ret;
            for (const auto& [qxi, ymap] : m_xmap & xi) {
                for (const auto& [qyi, qs] : ymap & yi) {
                    ret.push_back(region_t(rectangle_t(qxi, qyi), qs));
                }
            }
            return ret;
        }

      private:
        xmap_t m_xmap;
    };

    // A Rectangles which is indexed like a pixelized image or array.
    // The columns run along the X dimension, rows run along the Y.
    template<typename Value, typename Set = std::set<Value>>
    using Rectarray = Rectangles<Value, size_t, size_t, Set>;

    // Return a new Rectangles that is the result of mapping the input
    // Rectangles domain to a discretized "pixel" domain defined via
    // the binnings.  The bins of the binnings are interpreted as
    // pixel dimensions.  Any regions out of range of the Binnings are
    // truncated.
    template<typename Rec, typename Ret = Rectarray<typename Rec::value_t, typename Rec::set_t>>
    Ret pixelize(const Rec& inrec, const Binning xbins, const Binning& ybins)
    {
        using pxint_t = typename Ret::xinterval_t;
        using pyint_t = typename Ret::yinterval_t;
        Ret ret;

        const auto nxbins = xbins.nbins();
        const auto nybins = ybins.nbins();

        for (const auto& [xi, ymap] : inrec.xmap()) {
            auto pxi = pxint_t::right_open(std::clamp(xbins.bin(lower(xi)), 0, nxbins-1),
                                           std::clamp(xbins.bin(upper(xi)), 0, nxbins-1));
            
            for (const auto& [yi, val] : ymap) {
                auto pyi = pyint_t::right_open(std::clamp(ybins.bin(lower(yi)), 0, nybins-1),
                                               std::clamp(ybins.bin(upper(yi)), 0, nybins-1));
                ret.add(pxi, pyi, val);
            }
        }
        return ret;        
    }

    

}


#endif
