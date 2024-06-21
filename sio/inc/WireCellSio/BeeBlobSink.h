#ifndef WIRECELLSIO_BEEBLOBSINK
#define WIRECELLSIO_BEEBLOBSINK

#include "WireCellIface/IBlobSetSink.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/Bee.h"

namespace WireCell::Sio {

    /** BeeBlobSink sinks IBlobSets to Bee store
       
     */

    class BeeBlobSink
        : public Aux::Logger
        , public IConfigurable
        , public IBlobSetSink
        , public ITerminal
    {
    public:

        BeeBlobSink();
        virtual ~BeeBlobSink();

        virtual void finalize();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const IBlobSet::pointer& bs);
        
    private:

        /// Configure: outname
        std::string m_outname{"bee.zip"};

        /** Configuration: "geom"

            The canonical name of the Bee geometry.

            Required.

            Configuration: "type"

            The Bee algorithm type for display.

            Optional.
        */
        

        /** Configuration: "samplers"

            A string or array of string giving blob samplers.

            At least one required.
        */
        std::vector<IBlobSampler::pointer> m_samplers;

        /// Configure: run, sub, evt
        /// The evt will be added to the frame ident.
        int m_ident_offset{0};
        int m_last_ident{-1};

        Bee::Object m_obj;
        Bee::Sink m_sink;

        void flush(int ident = -1);

        size_t m_calls{0};

    };
}

#endif
