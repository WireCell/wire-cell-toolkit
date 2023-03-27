/** Convert from tensor data model to IFrame.
 */

#ifndef WIRECELL_AUX_TENSORFRAME
#define WIRECELL_AUX_TENSORFRAME

#include "WireCellIface/ITensorSetFrame.h"
#include "WireCellIface/IConfigurable.h"

#include <regex>

namespace WireCell::Aux {
    class TensorFrame : public WireCell::ITensorSetFrame,
                        public WireCell::IConfigurable
    {
      public:
        TensorFrame();
        virtual ~TensorFrame();
        
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
      private:

        /** Config: "dpre"

            The datapath at which the tensor representing a frame can
            be found.

            The string is interpreted as a regular expression (regex)
            and matched against the datapaths of the tensors.  This
            allows finding the frame tensor given a path that depends
            on the frame ident.
         */
        std::regex m_dpre;

        /** Config: decoding parameters.

            TensorFrame decodes its tensor input "enc" to produce the
            frame array "dec".

              dec = (enc - offset) / scale - baseline

            Defaults: offset=0, scale=1, baseline=0
        */
        double m_offset=0, m_scale=1, m_baseline=0;

        std::function<float(float)> m_transform;

    };
}

#endif
