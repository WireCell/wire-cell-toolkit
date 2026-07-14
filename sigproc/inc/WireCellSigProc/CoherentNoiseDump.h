#ifndef WIRECELLSIGPROC_COHERENTNOISEDUMP
#define WIRECELLSIGPROC_COHERENTNOISEDUMP

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/cnpy.h"

#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace WireCell {
    namespace SigProc {

        // Per-group capture for PDHD/PDVD coherent-noise sub validation.
        // All fields are detector-agnostic. Both PDHD::CoherentNoiseSub and
        // PDVD::CoherentNoiseSub populate one of these per apply() group when
        // debug_dump_path is non-empty; CoherentNoiseDumpWriter::write() emits
        // a single .npz per group.
        struct CoherentNoiseDump {
            // Group identity + geometry
            int apa = -1;
            int gid = -1;            // first-channel ident in the group
            int plane = -1;          // 0=U, 1=V, 2=W
            int nbin = 0;
            int res_offset = 0;
            std::vector<int> channels;          // per-channel ident, length nch

            // Knobs in scope (read from chndb at apply() entry)
            float decon_limit = 0.f;
            float decon_limit1 = 0.f;
            float roi_min_max_ratio = 0.f;
            float min_adc_limit = 0.f;
            float upper_adc_limit = 0.f;
            float upper_decon_limit = 0.f;
            float protection_factor = 0.f;
            int pad_front = 0;
            int pad_back = 0;
            std::string time_filter_name;
            std::string lf_tighter_filter_name;
            std::string lf_loose_filter_name;

            // Pre-SP capture (apply())
            WireCell::Waveform::realseq_t median;            // length nbin

            // Post-SP capture (filled by SignalProtection)
            WireCell::Waveform::realseq_t medians_decon_aligned;  // length nbin (zero if !decon_stage_ran)
            std::vector<uint8_t> signal_bool_raw;            // length nbin (threshold crossings, no pad)
            std::vector<uint8_t> signal_bool;                // length nbin (final, post-pad)
            float adc_threshold_chosen = 0.f;
            float decon_threshold_chosen = 0.f;
            float rms_adc = 0.f;
            float rms_decon = 0.f;
            float mean_adc = 0.f;
            float mean_decon = 0.f;
            uint8_t decon_stage_ran = 0;

            // Per-ROI (length nrois). roi_starts/ends are inclusive bin indices in
            // original-time. Padded extents are clamped to [0, nbin-1].
            std::vector<int32_t> roi_starts;
            std::vector<int32_t> roi_ends;

            // Per-ROI scalars derived from medians_decon (the median's acceptance
            // picture — what the user is tuning decon_limit1 / roi_min_max_ratio
            // against). Length nrois.
            std::vector<float> roi_max_median;
            std::vector<float> roi_min_median;
            std::vector<float> roi_ratio_median;            // |min|/max, or huge if max<=0
            std::vector<uint8_t> roi_accepted_median;       // (max>decon_limit1) && (|min|<max*ratio_cut)

            // Subtract_WScaling: per-channel scaling
            std::vector<float> scaling_coef;                 // length nch (parallel to channels)
            float ave_coef = 0.f;

            // Subtract_WScaling: per-(channel, ROI) decision matrix.
            // Stored row-major shape (nch, nrois). Same accept formula as
            // roi_accepted_median but applied per channel's own signal_roi_decon.
            std::vector<float> roi_max_per_ch;
            std::vector<float> roi_min_per_ch;
            std::vector<uint8_t> roi_accepted_per_ch;
        };

        struct CoherentNoiseDumpWriter {
            // Emit one .npz file under <root>/apa<apa>/<plane>_g<gid>.npz
            // containing every field of d. Caller has already gated on
            // m_dump_path being non-empty.
            static void write(const std::string& root, const CoherentNoiseDump& d)
            {
                static const char* plane_name[3] = {"U", "V", "W"};
                const std::string pname = (d.plane >= 0 && d.plane < 3) ? plane_name[d.plane] : "X";

                std::filesystem::path dir = std::filesystem::path(root) / ("apa" + std::to_string(d.apa));
                std::error_code ec;
                std::filesystem::create_directories(dir, ec);

                const std::string fpath = (dir / (pname + "_g" + std::to_string(d.gid) + ".npz")).string();
                // Truncate any prior file so each apply() writes a fresh archive.
                std::filesystem::remove(fpath, ec);

                auto save_vec = [&](const std::string& key, auto& v) {
                    using T = typename std::remove_reference_t<decltype(v)>::value_type;
                    if (v.empty()) {
                        T dummy{};
                        std::vector<size_t> shape = {0};
                        cnpy::npz_save(fpath, key, &dummy, shape, "a");
                    } else {
                        std::vector<size_t> shape = {v.size()};
                        cnpy::npz_save(fpath, key, v.data(), shape, "a");
                    }
                };
                auto save_scalar_f = [&](const std::string& key, float val) {
                    std::vector<size_t> shape = {1};
                    cnpy::npz_save(fpath, key, &val, shape, "a");
                };
                auto save_scalar_i = [&](const std::string& key, int32_t val) {
                    std::vector<size_t> shape = {1};
                    cnpy::npz_save(fpath, key, &val, shape, "a");
                };
                auto save_scalar_u8 = [&](const std::string& key, uint8_t val) {
                    std::vector<size_t> shape = {1};
                    cnpy::npz_save(fpath, key, &val, shape, "a");
                };
                auto save_string_as_chars = [&](const std::string& key, const std::string& s) {
                    if (s.empty()) {
                        char dummy = 0;
                        std::vector<size_t> shape = {0};
                        cnpy::npz_save(fpath, key, &dummy, shape, "a");
                        return;
                    }
                    std::vector<size_t> shape = {s.size()};
                    cnpy::npz_save(fpath, key, s.data(), shape, "a");
                };

                // Identity
                save_scalar_i("apa", d.apa);
                save_scalar_i("gid", d.gid);
                save_scalar_i("plane", d.plane);
                save_scalar_i("nbin", d.nbin);
                save_scalar_i("res_offset", d.res_offset);
                save_vec("channels", d.channels);

                // Knobs
                save_scalar_f("decon_limit", d.decon_limit);
                save_scalar_f("decon_limit1", d.decon_limit1);
                save_scalar_f("roi_min_max_ratio", d.roi_min_max_ratio);
                save_scalar_f("min_adc_limit", d.min_adc_limit);
                save_scalar_f("upper_adc_limit", d.upper_adc_limit);
                save_scalar_f("upper_decon_limit", d.upper_decon_limit);
                save_scalar_f("protection_factor", d.protection_factor);
                save_scalar_i("pad_front", d.pad_front);
                save_scalar_i("pad_back", d.pad_back);
                save_string_as_chars("time_filter_name", d.time_filter_name);
                save_string_as_chars("lf_tighter_filter_name", d.lf_tighter_filter_name);
                save_string_as_chars("lf_loose_filter_name", d.lf_loose_filter_name);

                // Waveforms
                save_vec("median", d.median);
                save_vec("medians_decon_aligned", d.medians_decon_aligned);
                save_vec("signal_bool_raw", d.signal_bool_raw);
                save_vec("signal_bool", d.signal_bool);

                // SP-stage scalars
                save_scalar_f("adc_threshold_chosen", d.adc_threshold_chosen);
                save_scalar_f("decon_threshold_chosen", d.decon_threshold_chosen);
                save_scalar_f("rms_adc", d.rms_adc);
                save_scalar_f("rms_decon", d.rms_decon);
                save_scalar_f("mean_adc", d.mean_adc);
                save_scalar_f("mean_decon", d.mean_decon);
                save_scalar_u8("decon_stage_ran", d.decon_stage_ran);

                // ROI
                save_vec("roi_starts", d.roi_starts);
                save_vec("roi_ends", d.roi_ends);
                save_vec("roi_max_median", d.roi_max_median);
                save_vec("roi_min_median", d.roi_min_median);
                save_vec("roi_ratio_median", d.roi_ratio_median);
                save_vec("roi_accepted_median", d.roi_accepted_median);

                // Subtract_WScaling
                save_vec("scaling_coef", d.scaling_coef);
                save_scalar_f("ave_coef", d.ave_coef);
                save_vec("roi_max_per_ch", d.roi_max_per_ch);
                save_vec("roi_min_per_ch", d.roi_min_per_ch);
                save_vec("roi_accepted_per_ch", d.roi_accepted_per_ch);
            }
        };

    }  // namespace SigProc
}  // namespace WireCell

#endif
