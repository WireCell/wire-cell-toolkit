/** A facade over the canonical optical "flash"/"light"/"flashlight" point clouds.
 *
 *  One optical flash: the singular flash quantities (time/value/ident/type) plus
 *  the per-optical-detector readouts (channel ident, time, value, error).  This is
 *  the generic, detector-agnostic read view shared by both the uboone path
 *  (root/UbooneClusterSource) and the SBND path (aux/FlashTensorToOpticalPCs) —
 *  both write the same schema.
 *
 *  Built by Grouping::flash_at()/flashes() (which own the flashlight-join walk);
 *  Cluster::get_flash() delegates here for a cluster's matched flash.  Read-only to
 *  callers; populated via friendship with Grouping.
 */

#ifndef WIRECELL_CLUS_FACADEFLASH
#define WIRECELL_CLUS_FACADEFLASH

#include <vector>
#include <cstddef>

namespace WireCell::Clus::Facade {

    class Grouping;

    class Flash {
        friend class Grouping;
        bool m_valid{false};
        double m_time{0}, m_value{0};
        int m_ident{-1}, m_type{-1};
        std::vector<int> m_idents;                       // per-OpDet channel id
        std::vector<double> m_times, m_values, m_errors; // per-OpDet readouts
        std::vector<int> m_cov_idents;                   // channels with coverage < 1
        std::vector<double> m_covs;                      // their covered fractions
    public:

        /// A "false" means there was no "flash" PC array and all values
        /// are invalid.  A "true" does not guarantee all values are valid.
        explicit operator bool() const { return m_valid; }

        /// Any "singular" methods are about the flash itself.

        /// Get the time of the flash.
        double time() const { return m_time; }

        /// Get the measurement of the flash
        double value() const { return m_value; }

        /// The ID of the flash
        int ident() const { return m_ident; }

        /// The type of the flash.
        int type() const { return m_type; }

        /// Any "plural" methods return per-optical-detector quantities.
        /// They will be empty() if the "light" and "flashlight" arrays are
        /// missing.  These vectors have the same size.

        // keep these return-by-value.

        /// Channel idents (optical-detector ids) of individual readouts.
        std::vector<int> idents() const { return m_idents; }
        /// Times of individual optical detector readouts.
        std::vector<double> times() const { return m_times; }
        /// Measurement values from optical detectors.
        std::vector<double> values() const { return m_values; }
        /// Measurement uncertainty from optical detectors.
        std::vector<double> errors() const { return m_errors; }

        /// Sparse per-flash readout-coverage rows (from the "flashcov" local
        /// PC written by Aux::FlashTensorToOpticalPCs when the light chain
        /// emitted coverage): channel ids with covered fraction < 1 of this
        /// flash's window, and those fractions.  Empty when the PC is absent
        /// (legacy archives) => everything fully covered.
        std::vector<int> cov_idents() const { return m_cov_idents; }
        std::vector<double> covs() const { return m_covs; }

        /// Per-channel measurement values indexed by OpDet id, zero-filled to
        /// `nchan` (a convenience for consumers needing a dense PE-by-channel
        /// vector, e.g. charge-light matching).  Out-of-range idents are skipped.
        std::vector<double> pes(int nchan) const {
            std::vector<double> v(nchan, 0.0);
            for (size_t i = 0; i < m_idents.size(); ++i) {
                const int ch = m_idents[i];
                if (ch >= 0 && ch < nchan) v[ch] = m_values[i];
            }
            return v;
        }
    };

}  // namespace WireCell::Clus::Facade

#endif
