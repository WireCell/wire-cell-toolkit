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

        /// A "false" means nothing was resolved, and every value below is the
        /// invalid default.  There are two producers of a Flash and each has its
        /// own reasons:
        ///
        ///   Grouping::flash_at(index)  -- no "flash" PC array, a negative
        ///     index, or an index past the last row.
        ///   Grouping::flash_by_gid(gid) -- no "opflash" PC, a negative gid, no
        ///     row carrying that gid, or a gid that names more than one flash
        ///     (see the uniqueness precondition on that method).
        ///
        /// A "true" means the singular values (time/value/ident/type) were read
        /// from real data.  It does NOT promise the plural per-OpDet vectors are
        /// filled (they are empty without the "light" and "flashlight" PCs), and
        /// for Cluster::get_flash() it does NOT promise the row is the flash
        /// that cluster actually matched -- see the caveat on
        /// Cluster::get_flash(), and prefer Cluster::get_matched_flash().
        ///
        /// The two producers also differ in what the fields MEAN: on the gid
        /// path ident() is the globally-unique gid rather than a per-input row
        /// id, type() is 0, errors() are 0 and cov_idents()/covs() are empty,
        /// because the "opflash" PC carries no type, saturation or coverage
        /// column.  time() and value() are identical on both.
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
