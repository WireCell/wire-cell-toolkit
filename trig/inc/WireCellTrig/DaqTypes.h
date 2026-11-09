/** Mirrors of the DUNE DAQ trigger data formats.

    These follow DUNE-DAQ trgdataformats develop (TP version 2, TA version 2,
    TC version 3) field for field, as plain members instead of bit fields, so
    that algorithm cores ported from triggeralgs keep the DAQ units and can be
    checked against it.  They do not depend on DAQ code.

    Units: time_* are timing-system (DTS) ticks, 16 ns at 62.5 MHz.  samples_*
    count ADC samples.  Converting one to the other needs the sample period:
    32 ticks per 512 ns sample, 31.25 per 500 ns sample.  See
    trig/docs/timing.org.
 */

#ifndef WIRECELLTRIG_DAQTYPES
#define WIRECELLTRIG_DAQTYPES

#include <cstdint>
#include <limits>
#include <vector>

namespace WireCell::Trig::Daq {

    using timestamp_t = uint64_t;
    using timestamp_diff_t = int64_t;
    using channel_t = uint32_t;
    using channel_diff_t = int32_t;
    using detid_t = uint8_t;
    using version_t = uint8_t;

    constexpr timestamp_t invalid_timestamp = std::numeric_limits<timestamp_t>::max();
    constexpr detid_t invalid_detid = std::numeric_limits<detid_t>::max();
    constexpr detid_t whole_detector = invalid_detid - 1;
    constexpr channel_t invalid_tp_channel = 0xFFFFFF;  // TP channel is 24 bits
    constexpr channel_t invalid_channel = std::numeric_limits<channel_t>::max();  // TA

    /// trgdataformats TriggerPrimitive, version 2.
    struct TriggerPrimitive {
        version_t version{2};
        uint8_t flag{0};
        detid_t detid{invalid_detid};
        channel_t channel{invalid_tp_channel};
        uint16_t samples_over_threshold{0xFFFF};
        timestamp_t time_start{invalid_timestamp};
        uint16_t samples_to_peak{0xFFFF};
        uint32_t adc_integral{0};
        uint16_t adc_peak{0};
    };

    /// trgdataformats TriggerActivityData, version 2.
    struct TriggerActivityData {
        enum class Type { kUnknown = 0, kTPC = 1, kPDS = 2 };
        enum class Algorithm {
            kUnknown = 0,
            kSupernova = 1,
            kPrescale = 2,
            kADCSimpleWindow = 3,
            kHorizontalMuon = 4,
            kMichelElectron = 5,
            kDBSCAN = 6,
            kPlaneCoincidence = 7,
            kChannelDistance = 8,
            kBundle = 9,
            kChannelAdjacency = 10,
            kSWIFT = 11,
            kProtoDUNEBSMWindow = 12,
        };
        version_t version{2};
        timestamp_t time_start{invalid_timestamp};
        timestamp_t time_end{invalid_timestamp};
        timestamp_t time_peak{invalid_timestamp};
        timestamp_t time_activity{invalid_timestamp};
        channel_t channel_start{invalid_channel};
        channel_t channel_end{invalid_channel};
        channel_t channel_peak{invalid_channel};
        uint64_t adc_integral{0};
        uint16_t adc_peak{0};
        detid_t detid{invalid_detid};
        Type type{Type::kUnknown};
        Algorithm algorithm{Algorithm::kUnknown};
    };

    /// triggeralgs TriggerActivity: the data plus copies of its TPs.
    struct TriggerActivity : public TriggerActivityData {
        std::vector<TriggerPrimitive> inputs;
    };

    /// trgdataformats TriggerCandidateData, version 3.
    struct TriggerCandidateData {
        enum class Type {
            kUnknown = 0, kTiming = 1, kTPCLowE = 2, kSupernova = 3, kRandom = 4, kPrescale = 5,
            kADCSimpleWindow = 6, kHorizontalMuon = 7, kMichelElectron = 8, kPlaneCoincidence = 9,
            kDBSCAN = 10, kChannelDistance = 11, kBundle = 12, kCTBFakeTrigger = 13, kCTBBeam = 14,
            kCTBBeamChkvHL = 15, kCTBCustomD = 16, kCTBCustomE = 17, kCTBCustomF = 18, kCTBCustomG = 19,
            kCTBBeamChkvHLx = 20, kCTBBeamChkvHxL = 21, kCTBBeamChkvHxLx = 22, kNeutronSourceCalib = 23,
            kChannelAdjacency = 24, kCIBFakeTrigger = 25, kCIBLaserTriggerP1 = 26, kCIBLaserTriggerP2 = 27,
            kCIBLaserTriggerP3 = 28, kCTBOffSpillSnapshot = 29, kCTBOffSpillCosmicJura = 30,
            kCTBOffSpillCRTCosmic = 31, kCTBBeamSpillStart = 32, kCTBBeamSpillSnapshot = 33,
            kCTBCustomC = 34, kCTBCustomPulseTrain = 35, kDTSPulser = 36, kDTSCosmic = 37,
            kSSPLEDCalibration = 38, kProtoDUNEBSMWindow = 39,
        };
        enum class Algorithm {
            kUnknown = 0, kSupernova = 1, kHSIEventToTriggerCandidate = 2, kPrescale = 3,
            kADCSimpleWindow = 4, kHorizontalMuon = 5, kMichelElectron = 6, kPlaneCoincidence = 7,
            kCustom = 8, kDBSCAN = 9, kChannelDistance = 10, kBundle = 11, kChannelAdjacency = 12,
        };
        version_t version{3};
        timestamp_t time_start{invalid_timestamp};
        timestamp_t time_end{invalid_timestamp};
        timestamp_t time_candidate{invalid_timestamp};
        detid_t detid{invalid_detid};
        Type type{Type::kUnknown};
        Algorithm algorithm{Algorithm::kUnknown};
    };

    /// triggeralgs TriggerCandidate: the data plus copies of its TAs (without
    /// their TPs).
    struct TriggerCandidate : public TriggerCandidateData {
        std::vector<TriggerActivityData> inputs;
    };

}  // namespace WireCell::Trig::Daq

#endif
