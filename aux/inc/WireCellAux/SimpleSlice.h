#ifndef WIRECELLAUX_SIMPLESLICE
#define WIRECELLAUX_SIMPLESLICE

#include "WireCellIface/ISlice.h"

namespace WireCell::Aux {

    struct SimpleSlice : public ISlice {
        IFrame::pointer m_frame{nullptr};
        int m_ident{0};
        double m_start{0}, m_span{0};
        ISlice::map_t m_activity{};

        SimpleSlice(IFrame::pointer frame, int ident, double start, double span, const ISlice::map_t& activity)
            : m_frame(frame), m_ident(ident), m_start(start), m_span(span), m_activity(activity) { }

        // Construct with no activity, user must fill in after construction.
        SimpleSlice(IFrame::pointer frame, int ident, double start, double span)
            : m_frame(frame), m_ident(ident), m_start(start), m_span(span) { }

        virtual ~SimpleSlice() {}

        // ISlice interface
        IFrame::pointer frame() const { return m_frame; }
        int ident() const { return m_ident; }
        double start() const { return m_start; }
        double span() const { return m_span; }
        map_t activity() const { return m_activity; }
        map_t& activity() { return m_activity; }
    };
}

#endif
