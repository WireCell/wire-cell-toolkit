#ifndef WIRECELL_SPNG_TENSORSELECTOR
#define WIRECELL_SPNG_TENSORSELECTOR

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchTensor.h"

#include <regex>

namespace WireCell::SPNG {

    /// TensorSelector provides rule-based selection and renaming of
    /// TDM-compliant ITorchTensor objects.
    ///
    /// Rules are provided as configuration as described below.
    ///
    /// This class is intended to be used as a (mixin) base class for an SPNG
    /// data flow graph node.  It is an IConfigurable which requires this
    /// interface class to be listed in the subclass's `WIRECELL_FACTORY()` CPP
    /// macro call.  If the subclass itself is configurable, it must include the
    /// TensorSelector's default_configuration() and forward the configuration
    /// object to the TensorSelector's configure() method.
    class TensorSelector : public IConfigurable {
    public:

        /// Configuration:
        ///
        /// tensor_selection :: an array of tensor selection rule objects.
        ///
        /// A selection rule object has the following attributes:
        ///
        /// - accept :: a regex to match against a tensor datapath.
        /// - reject :: a regex to match against a tensor datapath.
        ///
        /// The datapath of a tensor that is subject to selection is tested
        /// against the array of selection rule objects in the order of the
        /// array.  When the tensor is explicitly accepted (accept regex matches
        /// datapath) or rejected (reject regex matches datapath), the rule
        /// testing ceases.
        ///
        /// Policy for a tensor that matches no rules is left up to the user.
        ///
        ///
        /// tensor_renaming :: an array of tensor renaming rule objects.
        ///
        /// A renaming object has the following attributes
        ///
        /// - match :: a regular expression pattern applied to a tensor datapath.
        /// - replace :: a regular expression format to form the renamed datapath.
        ///
        /// The datapath of a tensor that is subject to renaming is tested
        /// against the array of renaming rule objects.  The first one that
        /// matches is applied and the remaining rules are ignored.
        ///
        /// Policy for a tensor that matches no rules is left up to the user.
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

        /// Return a selection result for a tensor.
        ///
        /// A kNoMatch means no accept and not reject matched.  Interpreting
        /// this is application defined.
        enum class SelectionResult { kNoMatch=0, kAccept=1, kReject=2 };
        SelectionResult select_tensor(const ITorchTensor::pointer ten) const;


        /// Apply datapath renaming rules to tensor.  If a rule matches, return
        /// a new tensor with the new datapath.  If no rule matches, return
        /// nullptr.
        ITorchTensor::pointer rename_tensor(const ITorchTensor::pointer ten) const;

    private:

        struct SelectionRule {
            std::regex accept, reject;
        };
        struct RenamingRule {
            std::regex match;
            std::string replace;
        };
        
        std::vector<SelectionRule> m_selection_rules;
        std::vector<RenamingRule> m_renaming_rules;


    };

}
#endif
