/** Synchronize frame streams.

    Synchronization is by frame ident number.

    It is assumed that each input stream provides frames with
    monotonically increasing frame ident numbers.

    The output frame stream is ordered by ident number.

    At any given time, once all input streams have at least one frame,
    the input frame with the smallest ident number is forwarded to the
    output.

    When an EOS is observed on any input then all frames on other
    inputs are collected until EOS is observed and these frames are
    output in ident-order followed by an EOS.

 */

#ifndef WIRECEL_AUX_FRAMESYNC
#define WIRECEL_AUX_FRAMESYNC

#include "WireCellIface/IFrameMerge.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Aux {

    class FrameSync : Aux::Logger, public IFrameMerge {
        bool m_flushing{false};
      public:

        FrameSync(const size_t multiplicity = 2);
        virtual ~FrameSync();
        virtual std::vector<std::string> input_types();

        // iqs: vector of input queues, mutable
        // oqs: tuple of output queues
        virtual bool operator()(input_queues& iqs, output_queues& oqs);

      private:
        size_t m_multiplicity {2};

        // recursive processing of input queues
        void flush(input_queues& iqs, output_queues& oqs);

        // input buffer
        input_queues m_iqs;
    };
}

#endif
