#ifndef WIRECELL_CONFIGMANAGER
#define WIRECELL_CONFIGMANAGER

#include "WireCellUtil/Configuration.h"

namespace WireCell {

    /** Bundle up some policy for handling configuration.
     *
     * Overall configuration object is a list of configuration
     * objects for each configurable instance.
     *
     * Each instance configuration is a dict with keys:
     *
     * - type gives the registered class name
     * - name gives an instance name, ("" by default if omitted)
     * - data gives a type-specific Configuration dictionary for the instance.
     *
     */
    class ConfigManager {
        Configuration m_top;

       public:
        ConfigManager();
        ~ConfigManager();

        /// Extend current list of configuration objects with more.
        /// Pass by value (kept for ABI compatibility with pre-built binaries
        /// such as cvmfs libWireCellApps.so); the body consumes `more` via
        /// per-entry std::move into m_top, so new callers should pass an
        /// rvalue (e.g. std::move(cfg)) to skip the caller-side deep-copy.
        void extend(Configuration more);

        // Add a fully-built configurable configuration for an instance, return its index
        int add(Configuration& cfg);

        // Add a configurable configuration by parts, return its index
        int add(Configuration& data, const std::string& type, const std::string& name = "");

        /// Return top-level, aggregate configuration
        const Configuration& all() const { return m_top; }

        Configuration at(int index) const;

        /// Return index of configuration for given class and instance
        /// names.  If not found, returned index == -1;
        int index(const std::string& type, const std::string& name = "") const;

        /// Return the number of configuration objects.
        int size() const { return m_top.size(); }

        /// Remove configuration at given index and return it.
        Configuration pop(int ind);

        /// Drop the bulk "data" sub-tree from every top-level entry, retaining
        /// only "type" and "name".  Useful after every IConfigurable has been
        /// constructed and configured: the resident Json::Value tree (which
        /// can be multiple GB for ICARUS-scale configs) is no longer needed
        /// for downstream operation since finalize() only consults type/name.
        void clear_data();

        /// Return a list of all known configurables
        typedef std::pair<std::string, std::string> ClassInstance;
        std::vector<ClassInstance> configurables() const;

    };

}  // namespace WireCell

#endif
