#ifndef WIRECELL_CLUS_MULTIALGBLOBCLUSTERING
#define WIRECELL_CLUS_MULTIALGBLOBCLUSTERING

#include "WireCellAux/Logger.h"
#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellUtil/Bee.h"

namespace WireCell::Clus {

    class MultiAlgBlobClustering
        : public Aux::Logger, public ITensorSetFilter, public IConfigurable, public ITerminal
    {
       public:
        MultiAlgBlobClustering();
        virtual ~MultiAlgBlobClustering() = default;

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual void finalize();

       private:

        /** Config: bee_zip
         *
         * The path to a zip file to fill with Bee JSON
         */
        Bee::Sink m_sink;
        int m_last_ident{-1};
        Bee::Points m_bee_img, m_bee_ld;
        Bee::Patches m_bee_dead;
        void flush(int ident = -1);
        void flush(WireCell::Bee::Points& bpts, int ident);

        bool m_save_deadarea{false};


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

        /** Config: "perf"
         *
         * If true, emit time/memory performance measures.  Default is false.
         */
        bool m_perf{false};

        // configurable parameters for dead-live clustering
        int m_dead_live_overlap_offset{2};

        // 
        IAnodePlane::pointer m_anode;
    };
}  // namespace WireCell::Clus

#endif
