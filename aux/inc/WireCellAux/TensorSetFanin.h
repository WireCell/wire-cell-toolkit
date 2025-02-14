/** An N-to-1 fan-in of ITensorSets.

    This will concatenate the input tensors in input port-order and output a
    single tensor set.

    The input ITensorSets themselves pose a potential point of collision
    w.r.t. the "ident" and "metadata".  See the docstring on the configuration
    parameters for how this collision may be handled.

*/
#ifndef WIRECELLAUX_TENSORSETFANIN
#define WIRECELLAUX_TENSORSETFANIN

#include "WireCellAux/Logger.h"

#include "WireCellIface/ITensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"

#include <vector>
#include <string>

namespace WireCell::Aux {

    class TensorSetFanin : public Aux::Logger, public ITensorSetFanin, public IConfigurable {
    public:
        TensorSetFanin();
        virtual ~TensorSetFanin();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // IFanin 
        virtual std::vector<std::string> input_types();
        virtual bool operator()(const input_vector& invec, output_pointer& out);

    private:
        /// Configuration: "ident_port" (optional, default=0)
        ///
        /// The output "ident" will be taken from the input on this port.
        int m_ident_port{0};

        // FIXME: add some "strategy" option to determine how ITensorSet
        // metadata is merged.  JSON merge-patch, take-one and dict.update().
        // It may also require a "metadata_order" option.
        // 
        // For now, the OUTPUT METADATA IS EMPTY!

        /// Configuration: "tensor_order" (optional, default=iota(nports))
        ///
        /// Provide an array of input port numbers.  If not provided, the input
        /// port order is assumed.  The tensors are concatenated in this order.
        /// The concatenation will honor any duplicate or missing port numbers.
        /// Ports count from 0.
        ///
        std::vector<size_t> m_tensor_order = {0};

        /// Configuration: "multiplicity" (required)
        ///
        /// Number of input ports.
        int m_multiplicity = {-1};

    };

}

#endif
