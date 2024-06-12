/**
   A source of IBlobSet from a MicroBooNE Wire-Cell Prototype files in ROOT
   format providing Trun, TC (aka "live") or TDC (aka "dead") TTrees.
 */

#ifndef WIRECELLROOT_UBOONEBLOBSETSOURCE
#define WIRECELLROOT_UBOONEBLOBSETSOURCE

#include "WireCellIface/IBlobSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/Logger.h"

#include "TFile.h"
#include "TTree.h"

namespace WireCell::Root {

    // Internal "bag" to hold various TTree addresses.
    class UbooneBlobSourceTrees;

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

            A string or array of string giving name(s) of file(s) to read for
            input.  Any missing files will be skipped with a log at warn level.
        */
        std::vector<std::string> m_input;
        
        /** Configuration: kind

            A string giving the name for the TTree holding clusters.

            "TC" for "live" clusters (default).
            "TDC" for "dead" clusters".

            Note: either will also load other info from the Trun TTree.
        */
        std::string m_kind{"live"};

        /** Configuration: views

            A list of strings of plane letters describing the view.

            For live kind, "uvw" outputs 3-view blobs with all live channels
            while "uv" outputs 2-view blobs with live channels where "w" view is
            ignored.

            For dead kind, "uv" produces 2-view blobs with dead channels and "w"
            channels ignored.

            To replicate typical output from WCT a live views list would be

              views=["uvw","uv","vw","wu"]

            And to replicate typical output from WCT a dead views list would be:

              views=["uv","vw","wu"]

            The default value is empty which is interpreted to me the above typical settings.
        */
        // Encode a view as bits {u,v,w}.
        std::vector<char> m_views;

        /** Configuration: anode

            Name the IAnodePlane component describing microboone.
        */
        IAnodePlane::pointer m_anode;
        IAnodeFace::pointer m_iface;

        // for logging
        size_t m_calls{0};

        // Connection to ROOT
        std::unique_ptr<UbooneBlobSourceTrees> m_data;

        // Mark if we have gone past our EOS.
        bool m_done{false};

        // Try to go to next entry, return false if input is exhausted.
        bool next();

        bool in_views(int bind);
        void bodge_channels(std::shared_ptr<Aux::SimpleSlice> slice, const RayGrid::Blob& blob, int bind);

        IFrame::pointer gen_frame();
        IChannel::pointer get_channel(int chanid);
        IBlobSet::pointer load_live();
        IBlobSet::pointer load_dead();
        std::pair<int,int> make_strip(const std::vector<int>& wire_in_plane_indices);
        RayGrid::Blob make_blob(int blob_in_tree_index);

    };

}

#endif
