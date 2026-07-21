/** Port of the SBND FrameShift derivation to bare-ROOT reco1 reading.
 *
 * prototype: sbndcode/Timing/FrameShift/FrameShift_module.cc
 * (SBN DocDB 46654, author Lan Nguyen).  The art module consumes the
 * decoded PTB and SPECTDC timestamp products plus the raw DAQ header and
 * derives, per event, the absolute (PPS-frame, ns) timestamps of the CRT
 * T1 reset, beam-gate opening and event-trigger frames.  Downstream
 * (yuhw's modified larwirecell OpFlashSource) a scalar reduction
 * FrameShiftInfo::FrameApplyAtCaf() of these frames is written into the
 * opflash tensor-set metadata as "frame_apply_at_caf" (ns).
 *
 * That reduction only exists in a local sbnobj modification (absent from
 * all public sbnobj branches as of 2026-07); the candidate implemented
 * here is frame_default - frame_etrig, to be validated against the
 * in-time-flash beam window (+0.3..1.9 us) and confirmed with the
 * authors.  See root/docs/sbnd-reco1-source.md.
 *
 * Only SBNDReco1OpFlashSource uses this header.
 */
#ifndef WIRECELLROOT_SBNDRECO1FRAMESHIFT
#define WIRECELLROOT_SBNDRECO1FRAMESHIFT

#include "WireCellRoot/SBNDReco1Products.h"

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace WireCell {
    namespace Root {
        namespace SBNDReco1 {

            constexpr uint64_t kInvalidTimestamp = std::numeric_limits<uint64_t>::max();

            // Defaults transcribed from
            // sbndcode/Timing/FrameShift/frameshift_sbnd_data.fcl
            struct FrameShiftConfig {
                uint64_t raw_ts_correction{367000};         // RawTSCorrection, ns (NTB offset)
                uint64_t shift_tdc_rwm2ptb_gate{1738};      // ShiftTdcRwm2PtbGate, ns
                uint64_t shift_tdc_etrig2ptb_etrig{257};    // ShiftTdcEtrig2PtbEtrig, ns
                uint64_t shift_tdc_crtt12ptb_crtt1{295};    // ShiftTdcCrtt12PtbCrtt1, ns
                std::vector<int> ptb_etrig_hlts{1, 2, 3, 4, 5, 14, 15, 16, 17};
                std::vector<int> beam_etrig_hlt{1, 2, 16};
                std::vector<int> offbeam_etrig_hlt{3, 4, 17};
                std::vector<int> xmuon_etrig_hlt{5, 14, 15};
                int beam_crtt1_hlt{20};
                int offbeam_crtt1_hlt{21};
                int beam_gate_hlt{26};
                int offbeam_gate_hlt{27};
            };

            struct FrameShiftResult {
                bool valid{false};
                std::string status;  // human-readable failure/summary
                bool is_beam{false}, is_offbeam{false}, is_xmuon{false};
                uint64_t raw_ts{kInvalidTimestamp};
                int hlt_etrig{-1};
                // Absolute PPS-frame timestamps (ns), kInvalidTimestamp when absent.
                uint64_t frame_crtt1{kInvalidTimestamp};
                uint64_t frame_gate{kInvalidTimestamp};
                uint64_t frame_etrig{kInvalidTimestamp};
                uint64_t frame_default{kInvalidTimestamp};
            };

            // prototype FrameShift_module.cc FindClosest()
            inline uint64_t find_closest(const std::vector<uint64_t>& timestamps, uint64_t reference)
            {
                uint64_t min_diff = kInvalidTimestamp;
                uint64_t closest = kInvalidTimestamp;
                for (const uint64_t& ts : timestamps) {
                    const uint64_t diff = reference > ts ? reference - ts : ts - reference;
                    if (diff < min_diff) {
                        min_diff = diff;
                        closest = ts;
                    }
                }
                return closest;
            }

            /** Derive the per-event frame shifts.
             *
             * raw_header_ts: artdaq::detail::RawEventHeader::timestamp (ns,
             * uncorrected -- the raw_ts_correction is applied here).
             * ptbs: decoded raw::ptb::sbndptb products (ptbdecoder).
             * tdcs: decoded sbnd::timing::DAQTimestamp products (tdcdecoder).
             */
            inline FrameShiftResult compute_frame_shift(uint64_t raw_header_ts,
                                                        const std::vector<raw::ptb::sbndptb>& ptbs,
                                                        const std::vector<sbnd::timing::DAQTimestamp>& tdcs,
                                                        const FrameShiftConfig& cfg)
            {
                FrameShiftResult res;

                // prototype GetRawTimestamp(): raw_ts = header timestamp - RawTSCorrection
                res.raw_ts = raw_header_ts - cfg.raw_ts_correction;

                // prototype GetTDCTimestamps(): bucket per channel.
                // ch0: CRT T1, ch1: BES, ch2: RWM, ch4: event trigger.
                std::vector<uint64_t> tdc_ch0, tdc_ch1, tdc_ch2, tdc_ch4;
                for (const auto& tdc : tdcs) {
                    switch (tdc.fChannel) {
                    case 0: tdc_ch0.push_back(tdc.fTimestamp); break;
                    case 1: tdc_ch1.push_back(tdc.fTimestamp); break;
                    case 2: tdc_ch2.push_back(tdc.fTimestamp); break;
                    case 4: tdc_ch4.push_back(tdc.fTimestamp); break;
                    }
                }

                // prototype GetPTBTimestamps(): unmask HLT trigger words into
                // per-bit (code, timestamp) entries; PTB timestamps count 20 ns
                // ticks (see decoder module) so scale by 20.
                std::vector<int> hlt_code;
                std::vector<uint64_t> hlt_ts;
                for (const auto& ptb : ptbs) {
                    for (const auto& trig : ptb.fHLTriggers) {
                        const uint64_t ts = trig.timestamp * uint64_t(20);
                        uint64_t val = trig.trigger_word;
                        for (int b = 0; b < 32; ++b) {
                            if (val & 0x1) {
                                hlt_code.push_back(b);
                                hlt_ts.push_back(ts);
                            }
                            val >>= 1;
                        }
                    }
                }
                if (hlt_code.empty()) {
                    res.status = "no PTB HLT timestamps";
                    return res;
                }

                // prototype FindETRIGs()
                uint64_t tdc_etrig_ts = kInvalidTimestamp;
                if (tdc_ch4.empty()) {
                    res.status = "no TDC ETRIG (ch4) timestamps";
                    return res;
                }
                tdc_etrig_ts = find_closest(tdc_ch4, res.raw_ts);

                int hlt_etrig = -1;
                uint64_t hlt_etrig_ts = kInvalidTimestamp;
                {
                    uint64_t min_diff = kInvalidTimestamp;
                    for (size_t i = 0; i < hlt_code.size(); ++i) {
                        bool is_etrig = false;
                        for (int code : cfg.ptb_etrig_hlts) {
                            if (hlt_code[i] == code) { is_etrig = true; break; }
                        }
                        if (!is_etrig) continue;
                        const uint64_t diff = res.raw_ts > hlt_ts[i] ? res.raw_ts - hlt_ts[i]
                                                                     : hlt_ts[i] - res.raw_ts;
                        if (diff < min_diff) {
                            min_diff = diff;
                            hlt_etrig = hlt_code[i];
                            hlt_etrig_ts = hlt_ts[i];
                        }
                    }
                }
                if (hlt_etrig < 0) {
                    res.status = "no HLT ETRIG timestamps";
                    return res;
                }
                res.hlt_etrig = hlt_etrig;

                // prototype DecideGlobalEtrigTimestamp(): TDC > PTB > raw header
                uint64_t global_etrig_ts = res.raw_ts;
                if (tdc_etrig_ts != kInvalidTimestamp) global_etrig_ts = tdc_etrig_ts;
                else if (hlt_etrig_ts != kInvalidTimestamp) global_etrig_ts = hlt_etrig_ts;

                // prototype DecideRelevantTDCTimestamps()
                uint64_t tdc_crtt1_ts = tdc_ch0.empty() ? kInvalidTimestamp : find_closest(tdc_ch0, global_etrig_ts);
                uint64_t tdc_rwm_ts = tdc_ch2.empty() ? kInvalidTimestamp : find_closest(tdc_ch2, global_etrig_ts);

                // prototype DecideRelevantPTBTimestamps(): stream from ETRIG HLT code
                int hlt_gate = -1, hlt_crtt1 = -1;
                for (int code : cfg.beam_etrig_hlt) {
                    if (hlt_etrig == code) {
                        hlt_gate = cfg.beam_gate_hlt;
                        hlt_crtt1 = cfg.beam_crtt1_hlt;
                        res.is_beam = true;
                        break;
                    }
                }
                if (!res.is_beam) {
                    for (int code : cfg.offbeam_etrig_hlt) {
                        if (hlt_etrig == code) {
                            hlt_gate = cfg.offbeam_gate_hlt;
                            hlt_crtt1 = cfg.offbeam_crtt1_hlt;
                            res.is_offbeam = true;
                            break;
                        }
                    }
                }
                if (!res.is_beam && !res.is_offbeam) {
                    for (int code : cfg.xmuon_etrig_hlt) {
                        if (hlt_etrig == code) {
                            res.is_xmuon = true;
                            break;
                        }
                    }
                }
                if (!res.is_beam && !res.is_offbeam && !res.is_xmuon) {
                    res.status = "ETRIG HLT " + std::to_string(hlt_etrig) + " matches no known stream";
                    return res;
                }

                uint64_t hlt_gate_ts = kInvalidTimestamp, hlt_crtt1_ts = kInvalidTimestamp;
                for (size_t i = 0; i < hlt_code.size(); ++i) {
                    if (hlt_code[i] == hlt_gate) hlt_gate_ts = hlt_ts[i];
                    if (hlt_code[i] == hlt_crtt1) hlt_crtt1_ts = hlt_ts[i];
                }

                // prototype produce() frame computation.
                if (res.is_beam || res.is_offbeam) {
                    // Frame CRT T1
                    if (tdc_crtt1_ts != kInvalidTimestamp)
                        res.frame_crtt1 = tdc_crtt1_ts - cfg.shift_tdc_crtt12ptb_crtt1;
                    else if (hlt_crtt1_ts != kInvalidTimestamp)
                        res.frame_crtt1 = hlt_crtt1_ts;
                    else
                        res.frame_crtt1 = 0;

                    // Frame beam gate: TDC RWM only for the beam stream
                    if (res.is_beam && tdc_rwm_ts != kInvalidTimestamp)
                        res.frame_gate = tdc_rwm_ts - cfg.shift_tdc_rwm2ptb_gate;
                    else if (hlt_gate_ts != kInvalidTimestamp)
                        res.frame_gate = hlt_gate_ts;
                    else
                        res.frame_gate = 0;
                }

                // Frame ETRIG
                if (tdc_etrig_ts != kInvalidTimestamp)
                    res.frame_etrig = tdc_etrig_ts - cfg.shift_tdc_etrig2ptb_etrig;
                else if (hlt_etrig_ts != kInvalidTimestamp)
                    res.frame_etrig = hlt_etrig_ts;
                else
                    res.frame_etrig = 0;

                // Default frame per stream: gate for (off)beam, ETRIG for xmuon.
                if (res.is_beam || res.is_offbeam)
                    res.frame_default = res.frame_gate;
                else
                    res.frame_default = res.frame_etrig;

                res.valid = true;
                res.status = res.is_beam ? "beam" : (res.is_offbeam ? "offbeam" : "xmuon");
                return res;
            }

        }  // namespace SBNDReco1
    }  // namespace Root
}  // namespace WireCell

#endif
