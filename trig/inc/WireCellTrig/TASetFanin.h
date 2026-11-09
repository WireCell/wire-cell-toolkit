/** Merge TA sets from several ports, such as several anodes, into one.

    This is a synchronous fanin: each call takes one set from every port.
    The sets must have the same coverage, as when the TA finders of all
    anodes run on frames of the same readout windows.  A mismatch is an
    error.  (Streams of different coverage would need an asynchronous merge
    with the earliest coverage end as watermark.)

    The merged TAs are sorted by time_start.  The output detector is the
    common one of the inputs, or -1 if they differ.  End of stream on any
    port ends the output stream.

    Configuration:

    - multiplicity :: the number of input ports (default 2).
 */

#ifndef WIRECELLTRIG_TASETFANIN
#define WIRECELLTRIG_TASETFANIN

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITriggerActivitySetFanin.h"

namespace WireCell::Trig {

    class TASetFanin : public Aux::Logger, public ITriggerActivitySetFanin, public IConfigurable {
       public:
        TASetFanin(size_t multiplicity = 2);
        virtual ~TASetFanin();

        virtual std::vector<std::string> input_types();
        virtual bool operator()(const input_vector& invec, output_pointer& out);

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

       private:
        size_t m_multiplicity;
        int m_count{0};
    };

}  // namespace WireCell::Trig

#endif
