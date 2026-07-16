#ifndef WIRECELL_MATCH_OPFLASH
#define WIRECELL_MATCH_OPFLASH

#include "WireCellClus/Facade_Flash.h"

#include <memory>
#include <set>
#include <vector>

namespace WireCell::Match {

    /// Per-channel PE-error model: PE_err = (PE < knee) ? floor : frac*PE.
    /// Defaults reproduce the historical hard-coded 0.3 rule; QLMatching passes
    /// its configured values so the model is tunable without a rebuild.
    struct PEErr {
        double floor = 0.3;  // PE_err for sub-knee channels
        double frac  = 0.3;  // fractional PE_err for channels at/above knee
        double knee  = 1.0;  // PE level (in PE) separating the two regimes
    };

    /// The matcher's per-flash working object: a thin adapter over the canonical
    /// optical flash (Clus::Facade::Flash) that adds the matching-specific
    /// conventions — synthesized per-channel PE_err (the 0.3 rule), the
    /// fired-channel list, and the time-ordering comparator used to keep the
    /// LASSO solve deterministic. It holds no tensor/PC knowledge of its own.
    class Opflash {
    public:
        using pointer = std::shared_ptr<Opflash>;

        /// Construct from the canonical flash facade: pulls time + the dense
        /// per-channel PE vector (pes(nchan)) and the flash ident, then
        /// synthesizes PE_err/fired via the matching convention below.
        /// pe_scale (optional, length nchan) multiplies the measured PE per
        /// channel before PE_err/total/fired are derived — a per-channel gain
        /// correction; nullptr => no scaling (byte-identical).
        Opflash(const WireCell::Clus::Facade::Flash& flash, double threshold, int nchan,
                const PEErr& pe_err = {}, const std::vector<double>* pe_scale = nullptr);

        /// Construct from a flash time and a per-channel PE vector
        /// (resized/zero-filled to nchan). The facade ctor delegates to this.
        /// PE_err is synthesized here (the PEErr rule), keeping that convention
        /// in one place. pe_scale as above.
        Opflash(double time, std::vector<double> pe, double threshold, int nchan,
                const PEErr& pe_err = {}, const std::vector<double>* pe_scale = nullptr);

        ~Opflash();

        void set_flash_id(int v) { flash_id = v; }
        void set_flash_type(int v) { type = v; }

        int    get_flash_id()  const { return flash_id; }
        double get_time()      const { return time; }
        double get_total_PE()  const { return total_PE; }
        const std::vector<double>& get_PEs() const { return PE; }
        double get_PE(int ch)     const { return PE[ch]; }
        double get_PE_err(int ch) const { return PE_err[ch]; }
        bool   get_fired(int ch)  const;
        int    get_num_fired()    const { return fired_channels.size(); }
        int    get_type()         const { return type; }
        double get_low_time()     const { return low_time; }
        double get_high_time()    const { return high_time; }
        int    get_num_channels() const { return m_nchan; }
        double get_threshold()    const { return m_threshold; }
        /// Per-flash per-channel saturation flag (DAPHNE rail overlap),
        /// carried on the light-PC "error" field by FlashTensorToOpticalPCs.
        /// All-zero unless the light chain ran with OpHitFinder
        /// flag_saturation; consumed by QLMatching use_saturation_flag.
        bool   get_sat(int ch)    const { return ch >= 0 && ch < (int)sat.size() && sat[ch]; }
        /// Per-flash per-channel readout-coverage fraction (self-trigger
        /// snippet livetime over this flash's window), carried on the
        /// sparse "flashcov" PC by FlashTensorToOpticalPCs.  1.0 when the
        /// PC is absent (legacy archives / full-stream channels); consumed
        /// by QLMatching use_coverage_flag.
        double get_cov(int ch)    const {
            if (cov.empty() || ch < 0 || ch >= (int)cov.size()) return 1.0;
            return cov[ch];
        }

        /// Widen PE_err to at least `err` on channels this flash never read
        /// out (get_cov < cov_min).  Such a channel measured "PE < the
        /// self-trigger threshold" (~1 PE on PDVD DAPHNE), not "PE == 0" to
        /// the tight PEErr floor.  Call after construction (the ctor fills
        /// cov after init() synthesizes PE_err).  Reaches the LASSO always;
        /// reaches the bundle chi2 only when pe_err_on_pred is off (else
        /// TimingTPCBundle derives the error from the predicted pe instead).
        /// No-op when err <= 0 or the light chain carried no coverage =>
        /// bit-identical by default.
        void inflate_nodata_err(double err, double cov_min);

    private:
        // Shared ctor body: fills PE/PE_err/total_PE/fired from a per-channel
        // PE vector (resized to nchan). flash_id is left 0 for callers to set.
        void init(double time, std::vector<double> pe, double threshold, int nchan,
                  const PEErr& pe_err, const std::vector<double>* pe_scale);

    protected:
        int    m_nchan;
        double m_threshold;

        int    type;
        int    flash_id;
        double low_time;
        double high_time;
        double time;
        double total_PE;

        std::vector<int>    fired_channels;
        std::vector<double> PE;
        std::vector<double> PE_err;
        std::vector<unsigned char> sat;  // per-channel saturation flags (empty = none)
        std::vector<float> cov;          // per-channel coverage fractions (empty = all 1)
    };

    struct OpFlashCompare {
        bool operator()(Opflash* a, Opflash* b) const {
            return a->get_time() < b->get_time();
        }
    };

    using OpflashSelection = std::vector<Opflash*>;
    using OpFlashSet       = std::set<Opflash*, OpFlashCompare>;

} // namespace WireCell::Match

#endif
