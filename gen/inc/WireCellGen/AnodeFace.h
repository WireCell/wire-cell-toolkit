#ifndef WIRECELLGEN_ANODEFACE
#define WIRECELLGEN_ANODEFACE

#include "WireCellIface/IAnodeFace.h"

namespace WireCell {
    namespace Gen {

        class AnodeFace : public IAnodeFace {
           public:
            AnodeFace(int ident, IWirePlane::vector planes, const BoundingBox& sensvol, int which, int aid);
            virtual ~AnodeFace();

            /// Return the ident number of this face.
            virtual int ident() const { return m_ident; }

            virtual int which() const { return m_which; }

            virtual int anode() const { return m_aid; }

            /// Return the number of wire planes in the given side
            virtual int nplanes() const { return m_planes.size(); }

            /// Return the wire plane with the given ident or nullptr if unknown.
            virtual IWirePlane::pointer plane(int ident) const;
            virtual IWirePlane::vector planes() const { return m_planes; }

            virtual BoundingBox sensitive() const { return m_bb; }

            virtual const RayGrid::Coordinates& raygrid() const { return m_coords; }

           private:
            int m_ident;
            IWirePlane::vector m_planes;
            BoundingBox m_bb;
            RayGrid::Coordinates m_coords;
            int m_which, m_aid;
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
