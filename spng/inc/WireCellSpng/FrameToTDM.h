/// Convert IFrame to ITorchTensorSet that is SPNG TDM compliant.

#ifndef WIRECELL_SPNG_FRAMETOTDM
#define WIRECELL_SPNG_FRAMETOTDM

#include "WireCellSpng/IFrameToTorchSet.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Waveform.h" // for ChannelMasks

namespace WireCell::SPNG {

    /// A FrameToTDM converts an IFrame to an ITorchTensorSet following SPNG
    /// torch data model.
    ///
    /// The IFrame and TDM "frame" data models tensors produced by FrameToTDM
    /// are not exactly equivalent for two reasons.
    ///
    /// Sparse vs dense: IFrame is, in general, a sparse set of ITrace.  ITrace may
    /// overlap, even when they are in the same tagged trace set.  How their
    /// overlap is accumulated is not specified by IFrame.  FrameToTDM will use
    /// additive accumulation to resolve any overlaps.
    ///
    /// Channel groups: For each tagged traces tag, FrameToTDM allows a number
    /// of "channel groups" to be defined.  The intersection of tagged traces
    /// and channels in a group determine the rows of a traces tensor and
    /// elements of a summaries tensor.  It is possible to specify group
    /// channels not represented in a tagged traces collection and vice versa.
    ///
    class FrameToTDM : public ContextBase, public Aux::Logger, public IFrameToTorchSet {
    public:

        FrameToTDM();
        virtual ~FrameToTDM();

        // IFunction
        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        

    private:

        /// Torch context configuration. 
        ///
        /// This is a ContextBase.  See that class for configuration.

        /// DATAPATHS:
        ///
        /// The FrameToTDM accepts user-defined datapaths assigned to each output tensor.
        ///
        /// A penultimate datapath for an output tensor is constructed as:
        ///
        ///   datapath = basepath + relpath
        ///
        /// Just prior to use, each datapath is formatted to interpolate any
        /// format codes like "{parameter}" with the following set of parameters:
        ///
        /// All tensors:
        ///
        /// - ident :: IFrame::ident() value
        ///
        /// In addition, channel masks take:
        ///
        /// - label :: the relevant label for a channel mask
        ///
        /// In addition, tenors made by "rules" take:
        ///
        /// - tag :: the rule's tag
        /// - index :: the 0-based index of the group in the group array.
        /// - part :: the relevant kind of frame "part" ("traces", "summaries", "chids").

        /// Configuration: "anode"
        ///
        /// The type:name of an anode plane component.
        IAnodePlane::pointer m_anode;

        /// Configuration: "basepath"
        ///
        /// The path APPENDED all individual relpaths with no special
        /// intervening character added.  If using filesystem-like paths, you
        /// must explicitly include a "/" either at end of your basepath or
        /// begin of all relpaths.
        std::string m_basepath = "/frames/{ident}/";

        /// Configuration: "frame_relpath"
        ///
        /// The path under basepath that the parent "frame" tensor is placed.
        std::string m_frame_relpath = "frame";

        /// Configuration: "chmasks"
        ///
        /// A map between a channel mask label ("bad", "noisy", etc) and a
        /// relative path pattern.  See DATAPATHS.  Only channel mask labels
        /// given are converted to tensors.
        std::map<std::string, std::string> m_chmasks;

        /// Configuration: "rules"
        ///
        /// An array of rule objects, each describing how derive tensors from the frame.
        ///
        /// Each rule object has the following attributes
        /// 
        /// - tag :: The tag indicating which "tagged traces" the rule consumes.
        ///   Multiple rules may have the same tag, which will lead to waveform
        ///   duplication.  The special tag "" (empty string) will instead apply
        ///   to all non-tagged traces.  Such an all-trace rule will not produce
        ///   any summaries tensor.
        ///
        /// - groups :: An array of "group objects" that defines how traces are
        ///   organized into "traces" tensors.  A group object has the following
        ///   attributes:
        ///
        ///   - wpids :: An array of WirePlaneId values (use wc.WirePlaneId() in
        ///     Jsonnet) that traces by their channels into a group.
        ///
        ///   - channels :: An array of channel IDs to select traces into a
        ///     group.  If both wpids and channels are given, a union of
        ///     channels are found.
        ///
        ///   - relpath :: A datapath relative to basepath under which
        ///     "traces/", "chids/", "summaries/" are placed.  This must be
        ///     unique across all group objects to avoid conflicts.  It is
        ///     recommended that the "tag" name be an element of this relative
        ///     path and a way to distinguish each group. 
        ///
        struct Group {
            std::string relpath; // relative to basepath
            std::string name;    // optional, "" by default
            std::set<int> wpids;
            std::set<int> channels;
        };
        struct Rule {
            std::string tag;    // traces tag
            std::vector<Group> groups;
        };

        std::vector<Rule> m_rules;
            
        // Map channel ID to a global ordering by major value wire-plane-id and
        // minor value wire-attachment-number.
        std::unordered_map<int, size_t> m_channel_order;

        // All the channels in a wireplane ID
        std::unordered_map<int, std::set<int> > m_wpid_channels;


        mutable size_t m_count{0};

    private:

        // Push comes to shove, I guess these can be specialized by a base
        // class.  For now, that's not the intention.

        /// Return the constituentless "frame" tensor
        ITorchTensor::pointer frame_itensor(const IFrame::pointer& iframe) const;

        /// Convert full channel mask map.
        ITorchTensor::vector chmask_tensors(const IFrame::pointer& iframe,
                                            const std::string& parent) const;
        torch::Tensor chmask_tensor(const Waveform::ChannelMasks& cms) const;

        
        /// Return the channel IDs that are in the intersection of have and group.
        ///
        /// The returned vector is sorted by increasing major value of wpid and
        /// minor value of wire-attachment-number (IChannel::index()).
        std::vector<int> group_channels(const std::vector<int>& have, const Group& group) const;

        /// Return trace indices in tagged set.  A tag of "" or "*" returns all
        /// indices simply enumerating those in the full traces vector.
        std::vector<size_t> tag_indices(const IFrame::pointer& iframe, const std::string& tag) const;

        /// Convert tensors according to the rules.
        ITorchTensor::vector rules_tensors(const IFrame::pointer& iframe, const std::string& parent) const;

    };
}

#endif
