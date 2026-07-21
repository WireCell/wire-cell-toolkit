#ifndef WIRECELL_SPNG_TENSORRENAMING
#define WIRECELL_SPNG_TENSORRENAMING

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchTensor.h"
#include "WireCellSpng/TensorIndex.h"

#include <regex>

namespace WireCell::SPNG {

    /// TensorRenaming provides rule-based renaming of TDM-compliant
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
    class TensorRenaming : public IConfigurable {
    public:

        /// The configuration object has a single attribute:
        ///
        /// - tensor_renaming :: an array of tensor renaming rule objects.
        ///
        /// A each renaming rule object in the tensor_renaming array has both of
        /// these attributes:
        ///
        /// - match :: a regular expression "pattern" string.
        ///
        /// - replace :: a regular expression "format" string.
        ///
        /// Each `match` regex from `tensor_renaming` array is applied in array
        /// order to the datapath of a tensor.  If successful the corresponding
        /// `replace` format is applied and any remainder of the
        /// `tensor_renaming` array is ignored.
        ///
        /// The user sets policy for when a tensor matches no rules.
        ///
        /// See the check-regex program provided by util/ as a way to validate
        /// your regex works are desired.  Some examples:
        ///
        /// @code{sh}
        /// $ check-regex "/path/to/dummy" '.*(dummy)' '/path/to/not-$1'
        /// /path/to/not-dummy
        ///
        /// $ check-regex "/path/to/dummy" '.*/dummy' '/path/to/smarty'
        /// /path/to/smarty
        ///
        /// $ check-regex "/path/to/dummy" 'dummy' 'smarty'
        /// /path/to/smarty
        ///
        /// $ check-regex "/path/to/dummy" '/(path)/(to)/(dummy)' '$3-$2-$1'
        /// dummy-to-path
        /// @endcode
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

        /// Return a tensor that has the datapath renaming rules applied.  If no
        /// rules apply, the original tensor is returned.  CAUTION: if this
        /// tensor is a parent of other tensors, this function will not update
        /// their "parent" metadata attribute.  See apply() for a way that
        /// applies consistent renaming across a set of tensors.
        ITorchTensor::pointer rename_tensor(const ITorchTensor::pointer ten) const;

        /// Return the name after passing through the renaming rules.  If no
        /// rules apply, the input datapath is returned.
        std::string rename(const std::string& datapath) const;

        /// Apply the renaming rules to the tensors in the index.  If a parent
        /// tensor has its datapath named, its children will have their "parent"
        /// attribute updated.
        TensorIndex apply(const TensorIndex& index) const;

    private:

        struct RenamingRule {
            std::regex match;
            std::string replace;
        };
        
        std::vector<RenamingRule> m_renaming_rules;


    };

}
#endif
