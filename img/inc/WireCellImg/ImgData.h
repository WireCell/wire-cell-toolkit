/**

   Some implementation of data interfaces used by INode components from wire-cell-img.

   Note: this header should NOT typically be #include'ed in other Wire
   Cell components which are not defined in wire-cell-img/src/.  It's
   exported as public to facilitiate writing unit tests.

*/

#include "WireCellIface/IFrame.h"
#include "WireCellIface/IWire.h"
#include "WireCellIface/IChannel.h"
#include "WireCellIface/IStripe.h"
#include "WireCellIface/IStripeSet.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/ISliceFrame.h"

namespace WireCell::Img::Data {

    class Slice : public ISlice {
        IFrame::pointer m_frame;
        map_t m_activity;
        int m_ident;
        double m_start, m_span;

      public:
        Slice(const IFrame::pointer& frame, int ident, double start, double span)
            : m_frame(frame)
            , m_ident(ident)
            , m_start(start)
            , m_span(span)
        {
        }

        virtual ~Slice();

        IFrame::pointer frame() const { return m_frame; }

        int ident() const { return m_ident; }
        double start() const { return m_start; }
        double span() const { return m_span; }
        map_t activity() const { return m_activity; }

        // These methods are not part of the ISlice interface and may be
        // used prior to interment in the ISlice::pointer.

        void sum(const IChannel::pointer& ch, value_t val) { m_activity[ch] += val; }
        void assign(const IChannel::pointer& ch, value_t val) { m_activity[ch] = val; }
        void merge(const ISlice::pointer& in) {m_activity.insert(in->activity().begin(),in->activity().end());}
    };

    // simple collection
    class SliceFrame : public ISliceFrame {
        ISlice::vector m_slices;
        int m_ident;
        double m_time;

      public:
        SliceFrame(const ISlice::vector& islices, int ident, double time)
            : m_slices(islices)
            , m_ident(ident)
            , m_time(time)
        {
        }
        virtual ~SliceFrame();

        int ident() const { return m_ident; }
        double time() const { return m_time; }
        ISlice::vector slices() const { return m_slices; }
    };

    class Stripe : public IStripe {
        int m_ident;
        vector_t m_values;

      public:
        Stripe(int ident)
            : m_ident(ident)
        {
        }
        virtual ~Stripe();

        int ident() const { return m_ident; }
        vector_t values() const { return m_values; }

        // these methods may be used prior to internment into IStripe::pointer

        void append(IChannel::pointer ich, value_t value) {
            m_values.push_back(make_pair(ich, value));
        }
    };
    class StripeSet : public IStripeSet {
        int m_ident;
        IStripe::vector m_stripes;

      public:
        StripeSet(int ident)
            : m_ident(ident)
        {
        }
        virtual ~StripeSet();

        int ident() const { return m_ident; }
        IStripe::vector stripes() const { return m_stripes; }

        // use before interning

        void push_back(const IStripe::pointer& s) { m_stripes.push_back(s); }
        size_t size() const { return m_stripes.size(); }
    };

}
