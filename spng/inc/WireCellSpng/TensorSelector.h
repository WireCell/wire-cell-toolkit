#ifndef WIRECELL_SPNG_TENSORSELECTOR
#define WIRECELL_SPNG_TENSORSELECTOR

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchTensor.h"
#include "WireCellSpng/TensorIndex.h"

#include <regex>

namespace WireCell::SPNG {

    /// TensorSelector provides rule-based selection of TDM-compliant
    /// ITorchTensor objects.
    ///
    /// Rules are provided as configuration as described below. 
    ///
    /// This class is intended to be used as a (mixin) base class or a data
    /// member for other component classes.  When used as a base class for a
    /// component, the IConfigurable interface class must be included in the
    /// subclass's `WIRECELL_FACTORY()` CPP macro call.  If the subclass is
    /// itself configurable it must marshal the configuration object from the
    /// base class default_configuration() and to the base class configure()
    /// methods.
    ///
    /// See also TensorSelector
    class TensorSelector : public IConfigurable {
    public:

        /// The configuration object has a single attribute:
        ///
        /// - tensor_selection :: an array of tensor selection rule objects.
        ///
        /// Each selection rule object in the `tensor_selection` array has one
        /// or both of the following attributes:
        ///
        /// - accept :: a regex to match against a tensor datapath.
        /// - reject :: a regex to match against a tensor datapath.
        ///
        /// A tensor's datapath is tested against the array of selection rule
        /// objects in the order of the array.  When the tensor is explicitly
        /// accepted (accept regex matches datapath) or rejected (reject regex
        /// matches datapath), the rule testing ceases.
        ///
        /// The user sets policy for the case that no rules match.
        ///
        /// See the check-regex program provided by util/ as a way to validate
        /// your regex works are desired.  Some examples:
        ///
        /// @code{sh}
        /// $ check-regex "/path/to/dummy" '.*/dummy' ; echo $status
        /// 0
        /// 
        /// $ check-regex "/path/to/dummy" '/path/to/dummy' ; echo $status
        /// 0
        ///
        /// $ check-regex "/path/to/dummy" 'nope' ; echo $status
        /// 1
        /// @endcode
        /// 
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

        /// Return a selection result for a tensor.
        ///
        /// A kNoMatch means no accept and not reject matched.  Interpreting
        /// this is application defined.
        enum class SelectionResult { kNoMatch=0, kAccept=1, kReject=2 };
        SelectionResult select_tensor(const ITorchTensor::pointer ten) const;

        /// Apply the selection rules to the parent tensors in the index
        /// (children of the root tree node).  Children live or die with their
        /// parent.  If keep_unselected is true (default) then any parent that
        /// is neither accepted nor rejected is considered accepted.  If false,
        /// it is considered rejected.
        TensorIndex apply(const TensorIndex& index, bool keep_unselected=true) const;



    private:

        struct SelectionRule {
            std::regex accept, reject;
        };

        std::vector<SelectionRule> m_selection_rules;

    };

}
#endif
