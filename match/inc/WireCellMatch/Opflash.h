#ifndef WIRECELL_MATCH_OPFLASH
#define WIRECELL_MATCH_OPFLASH

#include "WireCellIface/ITensorSet.h"

#include <memory>
#include <set>
#include <vector>

namespace WireCell::Match {

    class Opflash {
    public:
        using pointer = std::shared_ptr<Opflash>;

        /// Construct an Opflash from a 2D ITensor of doubles. ncol must be
        /// at least nchan+1. Column 0 is the flash time, columns 1..nchan
        /// hold per-channel PE.
        Opflash(const ITensor::pointer ten, int idx, double threshold, int nchan = 32);

        /// Construct an Opflash from a flash time and a per-channel PE vector
        /// (resized/zero-filled to nchan). This is the canonical-PC path; the
        /// tensor ctor above delegates to it. PE_err is synthesized here (the
        /// 0.3 rule), keeping that convention in one place.
        Opflash(double time, std::vector<double> pe, double threshold, int nchan);

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

    private:
        // Shared ctor body: fills PE/PE_err/total_PE/fired from a per-channel
        // PE vector (resized to nchan). flash_id is left 0 for callers to set.
        void init(double time, std::vector<double> pe, double threshold, int nchan);

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
