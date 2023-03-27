/*! A data interface to tensor data.
 *
 * This interface is sympathetic with Boost Multi Array.
 */

#ifndef WIRECELL_ITENSOR
#define WIRECELL_ITENSOR

#include "WireCellIface/IData.h"
#include "WireCellUtil/Configuration.h"
#include <vector>
#include <cstddef>  // std::byte
#include <typeinfo>
#include <typeindex>

namespace WireCell {

    class ITensor : public IData<ITensor> {
       public:
        /// Shape gives size of each dimension.  Size of shape give Ndim.
        typedef std::vector<size_t> shape_t;

        /// Storage order.  Empty implies C order.  If non-empty the
        /// vector holds the "majority" of the dimension.  C-order
        /// implies a vector of {1,0} which means if the array is
        /// accessed as array[a][b] "axis" 0 (indexed by "a") is the
        /// "major index" and "axis" 1 (indexed by "b") is the "minor
        /// index".  It is thus "row-major" ordering as the major
        /// index counts rows.  An array in fortran-order
        /// (column-major order) would be given as {0,1}.
        ///
        /// A note as this can be confusing: The "logical" rows and
        /// columns, eg when used in an Eigen array are independent
        /// from memory order.  An Eigen array is always indexed as
        /// arr(r,c).  Storage order only matters when, well, you
        /// access the array storage such as from Eigen array's
        /// .data() method - and indeed ITensor::data().
        typedef std::vector<size_t> order_t;

        /// The type of the element.
        virtual const std::type_info& element_type() const = 0;
        /// The size in bytes of an element
        virtual size_t element_size() const = 0;

        /// Helper to divine the dtype string, i4, f8, etc.
        virtual std::string dtype() const;

        /// The shape of the tensor.
        virtual shape_t shape() const = 0;

        /// The storage order of the tensor in the data array.
        virtual order_t order() const
        {
            return order_t();  // default is C storage order
        }

        /// The tensor store data array.
        virtual const std::byte* data() const = 0;
        /// The store size.  Size divided by product of shape gives
        /// word size
        virtual size_t size() const = 0;

        /// Optional metadata associated with the tensor
        virtual Configuration metadata() const { return Configuration(); }
    };
}  // namespace WireCell

#endif
