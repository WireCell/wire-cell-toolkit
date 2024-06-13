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

            A string "live" or "dead" describing what kind of blobs and slices to produce.

            Both require the Trun and TC TTrees.  In addition, "dead" requires
            the TDC TTree.  If the TTree named T_bad_ch exists it will also be
            loaded for both "live" and "dead".
        */
        std::string m_kind{"live"};

        /** Configuration: views

            A list of strings of plane letters describing the wire planes
            combinations to load.

            For the "live" kind, "uvw" outputs 3-view blobs with all live
            channels while "uv" outputs 2-view blobs with live channels where
            "w" view is ignored.  Etc for "vw" and "wu".  

            For dead kind, "uv" produces 2-view blobs with dead channels and "w"
            channels ignored.  Etc for "vw" and "wu".

            There is no ordering requirement on the letters in each string.

            To replicate typical output from WCT for "live" kind, the views list
            would be:

              views=["uvw","uv","vw","wu"]

            And to replicate typical output from WCT for "dead" kind, the views
            list would be:

              views=["uv","vw","wu"]

            The default value is an empty list which is interpreted to be the
            above full lists.
        */
        // Encode a view as OR'ed bits {u=1,v=2,w=4}
        std::vector<char> m_views;

        /** Configuration: anode

            Name the IAnodePlane component describing microboone.

            Required for looking up WCT channels given WCP wire indices (which
            for MB are identical to WCT wire-in-plane numbers).
        */
        IAnodePlane::pointer m_anode;
        IAnodeFace::pointer m_iface;

        // for logging
        size_t m_calls{0};

        // Our interface to ROOT
        std::unique_ptr<UbooneBlobSourceTrees> m_data;

        // Mark if we have gone past our EOS.
        bool m_done{false};

        // Try to go to next entry, return false if input is exhausted.
        bool next();

        // Return true if blob index is consistent with one of our configured "views".
        bool in_views(int bind);

        using SimpleSlicePtr = std::shared_ptr<Aux::SimpleSlice>;
        void bodge_channels(SimpleSlicePtr slice, const RayGrid::Blob& blob, int bind);

        IFrame::pointer gen_frame();
        IChannel::pointer get_channel(int chanid);

        // Map WCP timesliceId value to a slice
        using SliceMap = std::map<int, SimpleSlicePtr>;
        SliceMap load_slices();
        IBlobSet::pointer load_live();
        IBlobSet::pointer load_dead();
        std::pair<int,int> make_strip(const std::vector<int>& wire_in_plane_indices);
    };
}

#endif
