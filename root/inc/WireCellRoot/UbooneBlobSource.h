/**
   A source of IBlobSet from a MicroBooNE Wire-Cell Prototype files in ROOT
   format providing Trun, TC (aka "live"), TDC (aka "dead") and T_bad_ch TTrees.

   Each output blob set spans one slice (as it must).  Thus this source must be
   called multiple times to iterate through each ROOT TTree entry.

   Each blob set is constructed in exactly one seven supported cases:

   - 3-view live :: All three planes are categorized as *active*, none or *dummy* nor *masked*.
   - 2-view live :: Two planes are categorized as *active* and the third as *masked*, none are dummy.
   - 2-view dead :: No planes are *active*, two planes are *masked* and the third is *dummy*.

   See "views" and "kind" config parameters and the imaging-overview.org
   document section "Live vs dead" document for details of these categories.

 */

#ifndef WIRECELLROOT_UBOONEBLOBSETSOURCE
#define WIRECELLROOT_UBOONEBLOBSETSOURCE

#include "WireCellIface/IBlobSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/Logger.h"

#include "WireCellRoot/UbooneTTrees.h"

#include "TFile.h"
#include "TTree.h"

namespace WireCell::Root {

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

            Required.

            A string of two or three plane letters (any order of "u", "v" or
            "w") stating which planes are considered to have built the blobs.

            Blobs from file which are not consistent with kind+views are ignored.

            The "live" kind may take 3-views ("uvw") or any 2-views ("uv", "vw",
            "wu") and "dead" kind may take only 2-views.

            Default is "uvw" if kind is "live" else default is "uv".

        */
        // Encode a view as OR'ed bits {u=1,v=2,w=4}, can be in {3,5,6,7}.
        int m_views{0};
        // for dead, the third plane to get filled as dummy
        IWirePlane::pointer m_dummy{nullptr};
        std::vector<int> m_bodged; // for live/dead, plane indices that get bodged.

        /** Configuration: anode

            Name the IAnodePlane component describing microboone.

            Required for looking up WCT channels given WCP wire indices (which
            for MB are identical to WCT wire-in-plane numbers).
        */
        IAnodePlane::pointer m_anode;
        IAnodeFace::pointer m_iface;

        /** Configuration: frame_eos

            If true, emit an EOS after all blob sets for one ROOT TTree entry are output.

            This will require downstream to handle the stream restarting.

            Default is false.  Downstream may still identify an end of frame by
            the ident of the frame of the blobset stream changing.
        */
        bool m_frame_eos{false};

        // We don't even try to be just in time but for each TTree entry, fill a
        // queue of data to return over many calls.
        std::deque<IBlobSet::pointer> m_queue;

        // for logging
        size_t m_calls{0};

        // Our interface to the ROOT files and trees
        std::unique_ptr<UbooneTFiles> m_files;

        // Mark if we have gone past our EOS.
        bool m_done{false};

        // Return true if blob index is consistent with one of our configured "views".
        bool in_views(int bind);

        void bodge_activity(ISlice::map_t& activity, const RayGrid::Blob& blob);
        void dummy_activity(ISlice::map_t& activity);
        void fill_queue();

        IFrame::pointer gen_frame();
        IChannel::pointer get_channel(int chanid);

        void load_live();
        void load_dead();
        std::pair<int,int> make_strip(const std::vector<int>& wire_in_plane_indices);

        const double m_tick{0.5 * units::us};
        const ISlice::value_t m_bodge{0, 1e12};
    };
}

#endif
