#include "WireCellIface/ITrackSegmentSet.h"

namespace WireCell::Aux {

    class SimpleTrackSegmentSet : public ITrackSegmentSet {
        int m_ident;
        ITrackSegment::shared_vector m_segments;

       public:
        SimpleTrackSegmentSet(int ident, const ITrackSegment::vector& segments)
          : m_ident(ident)
          , m_segments(std::make_shared<ITrackSegment::vector>(segments.begin(), segments.end()))
        {
        }
        virtual ~SimpleTrackSegmentSet();
        virtual int ident() const { return m_ident; }
        virtual ITrackSegment::shared_vector segments() const { return m_segments; }
    };

}  // namespace WireCell::Aux
