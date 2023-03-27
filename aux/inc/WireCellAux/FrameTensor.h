#ifndef WIRECELL_AUX_FRAMETENSOR
#define WIRECELL_AUX_FRAMETENSOR

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameTensorSet.h"
#include "WireCellAux/TensorDM.h"

namespace WireCell::Aux {
    class FrameTensor : public WireCell::IConfigurable,
                        public WireCell::IFrameTensorSet
    {
      public:
        FrameTensor();
        virtual ~FrameTensor();
        
        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);
        
        /** Config: "datapath"

            Set the datapath for the tensor representing the frame.
            The string may provide a %d format code which will be
            interpolated with the frame's ident number.
         */
        std::string m_datapath = "frames/%d";

        /** Config: "mode".

            Set the mode to map frame to tensor to one of "sparse",
            "monolithic" or "tagged".
        */
        TensorDM::FrameTensorMode m_mode;

        /** Config: "digitize".

            If true, store trace as short ints, otherwise as floats.
        */
        bool m_digitize{false};

        /** Config: encoding parameters.

            FrameTensor encodes its frame input "dec" to produce
            tensor array "enc".

              enc = (dec + baseline) * scale + offset

            Defaults: offset=0, scale=1, baseline=0
        */
        double m_offset=0, m_scale=1, m_baseline=0;

        std::function<float(float)> m_transform;

    };
}

#endif
