#ifndef WIRECELL_CLUS_POINTTREEMERGING
#define WIRECELL_CLUS_POINTTREEMERGING


#include "WireCellAux/Logger.h"
#include "WireCellIface/ITensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITerminal.h"

namespace WireCell::Clus {

    class PointTreeMerging
        : public Aux::Logger, public ITensorSetFanin, public IConfigurable, public ITerminal
    {
       public:
        PointTreeMerging();
        virtual ~PointTreeMerging() = default;

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> input_types();

        // ITensorSetFanin, a PCTree for each input ITensorSet
        virtual bool operator()(const input_vector& invec, output_pointer& out);

        virtual void finalize();

       private:

        /** Config: "inpath"
         *
         * The datapath for the input point graph data.  This may be a
         * regular expression which will be applied in a first-match
         * basis against the input tensor datapaths.  If the matched
         * tensor is a pcdataset it is interpreted as providing the
         * nodes dataset.  Otherwise the matched tensor must be a
         * pcgraph.
         */
        std::string m_inpath{".*"};

        /** Config: "outpath"
         *
         * The datapath for the resulting pcdataset.  A "%d" will be
         * interpolated with the ident number of the input tensor set.
         */
        std::string m_outpath{""};

        /** Config: "perf"
         *
         * If true, emit time/memory performance measures.  Default is false.
         */
        bool m_perf{true};

        // Count how many times we are called
        size_t m_count{0};
        size_t m_multiplicity{0};
    };
}

#endif
