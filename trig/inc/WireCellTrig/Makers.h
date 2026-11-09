/** Base classes for trigger activity (TA) and candidate (TC) makers.

    These reproduce the machinery of the DUNE-DAQ triggeralgs develop
    TriggerActivityMaker and TriggerCandidateMaker: preprocess (input filter),
    process (the algorithm), postprocess (prescale) and flush.  They work on
    the DAQ data mirrors of DaqTypes.h, in DAQ units, and use no WCT interface
    types.

    Differences from triggeralgs, all deliberate:

    - Prescale counts only the objects appended by the current call.
      triggeralgs re-scans the whole output vector, which miscounts when a
      caller accumulates outputs over calls (as dunetrigger does).

    - store_inputs=false empties the "inputs" of the emitted objects.  The
      algorithms still see their inputs internally.  triggeralgs always stores.

    - flush() is part of the contract: a stream driver calls it when it knows
      no input earlier than "until" will come, and at end of stream.  The DAQ
      (develop) never calls flush().
 */

#ifndef WIRECELLTRIG_MAKERS
#define WIRECELLTRIG_MAKERS

#include "WireCellTrig/DaqTypes.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace WireCell::Trig {

    struct MakerConfig {
        /// Keep one in "prescale" of the made objects.
        uint64_t prescale{1};
        /// Whether emitted objects keep copies of their inputs.
        bool store_inputs{true};
    };

    class TAMaker {
       public:
        using input_type = Daq::TriggerPrimitive;
        using output_type = Daq::TriggerActivity;
        using output_vector = std::vector<output_type>;

        struct Config : public MakerConfig {
            /// TPs longer than this are ignored (triggeralgs
            /// max_samples_over_threshold).
            uint32_t max_samples_over_threshold{std::numeric_limits<uint32_t>::max()};
        };

        virtual ~TAMaker();

        void set_common(const Config& cfg) { m_common = cfg; }
        const Config& common() const { return m_common; }

        /// Feed one TP, append any TAs made.
        void operator()(const input_type& tp, output_vector& out);

        /// Emit what can be emitted given no TP before "until" will come.
        void flush(Daq::timestamp_t until, output_vector& out);

       protected:
        virtual bool preprocess(const input_type& tp);
        virtual void process(const input_type& tp, output_vector& out) = 0;
        virtual void do_flush(Daq::timestamp_t until, output_vector& out);

       private:
        Config m_common;
        uint64_t m_count{0};
        void postprocess(output_vector& out, size_t first);
    };

    class TCMaker {
       public:
        using input_type = Daq::TriggerActivity;
        using output_type = Daq::TriggerCandidate;
        using output_vector = std::vector<output_type>;

        struct Config : public MakerConfig {
            /// The TC type to give to the made TCs (triggeralgs tc_type_name).
            Daq::TriggerCandidateData::Type tc_type{Daq::TriggerCandidateData::Type::kUnknown};
        };

        virtual ~TCMaker();

        void set_common(const Config& cfg) { m_common = cfg; }
        const Config& common() const { return m_common; }

        void operator()(const input_type& ta, output_vector& out);
        void flush(Daq::timestamp_t until, output_vector& out);

       protected:
        virtual bool preprocess(const input_type& ta);
        virtual void process(const input_type& ta, output_vector& out) = 0;
        virtual void do_flush(Daq::timestamp_t until, output_vector& out);

       private:
        Config m_common;
        uint64_t m_count{0};
        void postprocess(output_vector& out, size_t first);
    };

}  // namespace WireCell::Trig

#endif
