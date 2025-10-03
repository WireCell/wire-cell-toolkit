/// Convert IFrame to ITorchTensorSet that is SPNG TDM compliant.

#ifndef WIRECELL_SPNG_FRAMETOTDM
#define WIRECELL_SPNG_FRAMETOTDM

#include "WireCellSpng/IFrameToTorchSet.h"
#include "WireCellSpng/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Waveform.h" // for ChannelMasks

#include "boost/filesystem.hpp"

namespace WireCell::SPNG {

    /// A FrameToTdm converts an IFrame to an ITorchTensorSet following SPNG
    /// torch data model.
    ///
    /// The IFrame and Tdm "frame" data models tensors produced by FrameToTdm
    /// are not exactly equivalent for two reasons.
    ///
    /// Sparse vs dense: IFrame is, in general, a sparse set of ITrace.  ITrace may
    /// overlap, even when they are in the same tagged trace set.  How their
    /// overlap is accumulated is not specified by IFrame.  FrameToTdm will use
    /// additive accumulation to resolve any overlaps.
    ///
    /// Channel groups: For each tagged traces tag, FrameToTdm allows a number
    /// of "channel groups" to be defined.  The intersection of tagged traces
    /// and channels in a group determine the rows of a traces tensor and
    /// elements of a summaries tensor.  It is possible to specify group
    /// channels not represented in a tagged traces collection and vice versa.
    ///
    class FrameToTdm : public Logger,
                       public virtual IConfigurable,
                       public virtual IFrameToTorchSet {
    public:

        FrameToTdm();
        virtual ~FrameToTdm();

        // IFunction
        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
        /// CONFIGURATION:
        ///
        /// FrameToTdm has sophisticated configuration but most uses can rely on
        /// defaults for most parameters.  See the document
        /// spng/docs/frametotdm.org for details.

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        

    private:

        /// The basepath for the datapaths of all output tensors.
        boost::filesystem::path m_basepath = "/frames/{ident}";

        /// Configuration: "frame_relpath"
        ///
        /// The path under basepath to place the parent "frame" tensor.
        boost::filesystem::path m_frame_relpath = "frame";

        /// Configuration: "anode"
        ///
        /// The type:name of an anode plane component.
        IAnodePlane::pointer m_anode;

        /// Configuration: "chmasks"
        ///
        /// Each "label" key in the "channel mask map" is converted to a
        /// "chmasks" TDM tensor.  The datapath of each label may be customized
        /// by providing a "chmasks" object that maps label attribute key to a
        /// "relpath" string.  During conversion, the empty string ("") label
        /// entry will be used if the label is not found in the "chmasks"
        /// configuration object.  This default relpath is "chmasks/{label}".
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
        ///   any summaries tensor.  Note, when "tag" is used in formatting a
        ///   path, a tag of "" is spelled as "null".
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
        ///   - relpath :: A datapath relative to basepath under which the parts
        ///     "traces", "chids", "summaries" are placed.  This must be unique
        ///     across all group objects to avoid conflicts.  It is recommended
        ///     to include the '{tag}' if not empty and '{index}' if more than
        ///     one group is in the array and '{datatype}' if multiple parts are
        ///     extracted.
        ///
        struct Group {
            /// The path relative to basepath for each part tensor in the group.
            /// See spng/docs/frametotdm.org
            boost::filesystem::path relpath {
                "tags/{tag}/rules/{rule}/groups/{group}/{datatype}"
            };
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

        // Helper for some datatypes with common patterns.  All torch::Tensor go
        // through here.
        ITorchTensor::pointer make_datatype(const std::string& datatype,
                                            const boost::filesystem::path& relpath,
                                            torch::Tensor ten, Configuration md) const;

    };
}

#endif
