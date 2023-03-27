/** SimpleBlob
 */
#ifndef WIRECELLAUX_SIMPLEBLOB
#define WIRECELLAUX_SIMPLEBLOB

#include "WireCellIface/IBlob.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellUtil/RayGrid.h"

namespace WireCell::Aux {

    class SimpleBlob : public IBlob {
       public:
        SimpleBlob(int ident, float value, float uncertainty,
                   const RayGrid::Blob& shape, ISlice::pointer slice,
                   IAnodeFace::pointer face)
          : m_ident(ident)
          , m_value(value)
          , m_uncertainty(uncertainty)
          , m_shape(shape)
          , m_slice(slice)
          , m_face(face)
        {
        }
        virtual ~SimpleBlob();

        int ident() const { return m_ident; }

        float value() const { return m_value; }

        float uncertainty() const { return m_uncertainty; }

        IAnodeFace::pointer face() const { return m_face; }

        ISlice::pointer slice() const { return m_slice; }

        const RayGrid::Blob& shape() const { return m_shape; }

       private:
        int m_ident;
        float m_value;
        float m_uncertainty;
        RayGrid::Blob m_shape;
        ISlice::pointer m_slice;
        IAnodeFace::pointer m_face;
    };

    class SimpleBlobSet : public IBlobSet {
       public:
        SimpleBlobSet(int ident, const ISlice::pointer& slice,
                      const IBlob::vector& blobs = {})
          : m_ident(ident)
          , m_slice(slice)
          , m_blobs(blobs)
        {
        }
        virtual ~SimpleBlobSet();

        virtual int ident() const { return m_ident; }

        virtual ISlice::pointer slice() const { return m_slice; }

        virtual IBlob::vector blobs() const { return m_blobs; }

        int m_ident;
        ISlice::pointer m_slice;
        IBlob::vector m_blobs;
    };

}  // namespace WireCell::Aux

#endif

