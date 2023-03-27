#ifndef WIRECELLIMG_POINTGRAPHPROCESSOR
#define WIRECELLIMG_POINTGRAPHPROCESSOR

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IPointGraphVisitor.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Img {

    class PointGraphProcessor : public Aux::Logger, public ITensorSetFilter, public IConfigurable {
      public:
        PointGraphProcessor();
        virtual ~PointGraphProcessor() = default;

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;         

        virtual bool operator()(const input_pointer& in, output_pointer& out);
         
      private:
        // Count how many times we are called
        size_t m_count{0};

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

        /** Config: "visitors"
         *
         * An array of string giving the type/name of
         * IPointGraphVisitor instances.
         */
        std::vector<IPointGraphVisitor::pointer> m_visitors;
         
    };
}

#endif
