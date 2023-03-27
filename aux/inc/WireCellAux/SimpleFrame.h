#ifndef WIRECELLAUX_SIMPLEFRAME
#define WIRECELLAUX_SIMPLEFRAME

#include "WireCellIface/IFrame.h"
#include "WireCellUtil/Units.h"

#include <map>

namespace WireCell::Aux {

    /** A simple frame.
     *
     * This is is nothing more than a bag of data.
     */
    class SimpleFrame : public IFrame {
       public:
        SimpleFrame(int ident, double time, const ITrace::vector& traces, double tick = 0.5 * units::microsecond,
                    const Waveform::ChannelMaskMap& cmm = Waveform::ChannelMaskMap());

        SimpleFrame(int ident, double time, ITrace::shared_vector traces, double tick = 0.5 * units::microsecond,
                    const Waveform::ChannelMaskMap& cmm = Waveform::ChannelMaskMap());

        ~SimpleFrame();
        virtual int ident() const;
        virtual double time() const;
        virtual double tick() const;

        virtual ITrace::shared_vector traces() const;
        virtual Waveform::ChannelMaskMap masks() const;

        virtual const tag_list_t& frame_tags() const;
        virtual const tag_list_t& trace_tags() const;
        virtual const trace_list_t& tagged_traces(const tag_t& tag) const;
        virtual const trace_summary_t& trace_summary(const tag_t& tag) const;

        // Before this object is interned to a shared vector the
        // creator may call the following non-const methods:

        // Add a frame tag.
        void tag_frame(const tag_t& tag);

        // Tag a subset of traces with optional trace summary

        template<typename Indices, typename Summary>
        void tag_traces(const tag_t& tag, const Indices& indices, const Summary& summary)
        {
            m_trace_info[tag] = SimpleTraceInfo(indices, summary);
            m_trace_tags.clear();
            for (auto& it : m_trace_info) {
                m_trace_tags.push_back(it.first);
            }
        }

        template<typename Indices>
        void tag_traces(const tag_t& tag, const Indices& indices)
        {
            m_trace_info[tag] = SimpleTraceInfo(indices);
            m_trace_tags.clear();
            for (auto& it : m_trace_info) {
                m_trace_tags.push_back(it.first);
            }
        }

        // old form
        void tag_traces(const tag_t& tag, const IFrame::trace_list_t& indices,
                        const IFrame::trace_summary_t& summary = IFrame::trace_summary_t(0))
        {
            tag_traces<IFrame::trace_list_t,IFrame::trace_summary_t>(tag, indices, summary);
        }

       private:
        int m_ident;
        double m_time, m_tick;
        ITrace::shared_vector m_traces;
        Waveform::ChannelMaskMap m_cmm;

        IFrame::tag_list_t m_frame_tags, m_trace_tags;

        struct SimpleTraceInfo {
            IFrame::trace_list_t indices;
            IFrame::trace_summary_t summary;
            SimpleTraceInfo();

            template<typename Indices>
            explicit SimpleTraceInfo(const Indices& inds)
                : indices(inds.begin(), inds.end())
            {
            }
            template<typename Indices, typename Summary>
            SimpleTraceInfo(const Indices& inds, const Summary& sums)
                : indices(inds.begin(), inds.end())
                , summary(sums.begin(), sums.end())
            {
            }
        };
        std::map<IFrame::tag_t, SimpleTraceInfo> m_trace_info;

        const SimpleTraceInfo& get_trace_info(const IFrame::tag_t& tag) const;
    };

}  // namespace WireCell::Aux

#endif
