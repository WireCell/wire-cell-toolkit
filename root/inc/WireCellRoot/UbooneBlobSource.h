/**
   A source of IBlobSet from a MicroBooNE Wire-Cell Prototype file in ROOT
   format providing Trun, TC (or TDC) TTrees.

 */

#ifndef WIRECELLROOT_UBOONEBLOBSETSOURCE
#define WIRECELLROOT_UBOONEBLOBSETSOURCE

#include "WireCellIface/IBlobSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellAux/Logger.h"

#include "TFile.h"
#include "TTree.h"

namespace WireCell::Root {

    struct UbooneBlobSourceTrees;

    class UbooneBlobSource : public Aux::Logger,
                                    public IBlobSetSource,
                                    public IConfigurable
    {
    public:
        UbooneBlobSource();
        virtual ~UbooneBlobSource();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        virtual bool operator()(IBlobSet::pointer& blobset);
    private:

        /** Configuration: input

            A string or array of string giving name(s) of file(s)
            to read for input.
        */
        std::vector<std::string> m_input;
        
        /** Configuration: kind

            A string giving the name for the TTree holding clusters.

            "TC" for "live" clusters (default).
            "TDC" for "dead" clusters".
            Note: either will also load other info from the Trun TTree.

        */
        std::string m_kind{"live"};

        /** Configuration: anode

            Name the IAnodePlane component.
        */
        IAnodePlane::pointer m_anode;
        IAnodeFace::pointer m_iface; // first face

        size_t m_calls{0};

        Long64_t m_entry{0};
        std::unique_ptr<TFile> m_tfile;
        TTree* m_trun;
        TTree* m_tc;

        std::unique_ptr<UbooneBlobSourceTrees> m_data;

        bool next();
        IFrame::pointer gen_frame();
        IChannel::pointer get_channel(int chanid);
        IBlobSet::pointer load();

    };

}

#endif
