/** JSON configuration of the trig algorithm cores, for the WCT nodes.

    The cores themselves take plain structs.  Keys of the ported DAQ
    algorithms are those of triggeralgs, in its units (DTS ticks, ADC, cm).

    Note, WireCell::get<T>() silently returns the default for types other
    than bool, int, int64_t, float, double and string, so these read through
    those and cast.
 */

#ifndef WIRECELLTRIG_TRIGCONFIG
#define WIRECELLTRIG_TRIGCONFIG

#include "WireCellTrig/BurstMaker.h"
#include "WireCellTrig/Makers.h"
#include "WireCellTrig/SWIFTMaker.h"
#include "WireCellTrig/SupernovaMakers.h"
#include "WireCellTrig/TPGChannel.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell::Trig {

    TPGConfig tpg_config_from_json(const Configuration& cfg);
    Configuration tpg_config_to_json(const TPGConfig& tpg);

    /// Keys: prescale, max_samples_over_threshold.  (store_inputs is handled
    /// by the nodes.)
    TAMaker::Config tamaker_config_from_json(const Configuration& cfg);

    /// Keys: prescale, tc_type (an integer, DAQ TriggerCandidateData::Type).
    TCMaker::Config tcmaker_config_from_json(const Configuration& cfg);

    /// Keys: time_tolerance, channel_tolerance, both, ticks_per_sample.
    TAMakerSupernova::Config supernova_ta_config_from_json(const Configuration& cfg,
                                                           double default_ticks_per_sample = 32);

    /// Keys: time_window, threshold, hit_threshold.
    TCMakerSupernova::Config supernova_tc_config_from_json(const Configuration& cfg);

    /// Keys: window_length, inspect_energy_threshold_sadc,
    /// accept_energy_threshold_sadc, min_adc_peak, min_samples_over_threshold,
    /// cm_per_tick, wire_pitch, min_samples, epsilon, cluster_energy_cut_sadc.
    TAMakerSWIFT::Config swift_config_from_json(const Configuration& cfg);

    /// Keys, in WCT units: window, threshold, holdoff, tc_type, algorithm,
    /// and the selection: algorithms (array of int), min_inputs, max_inputs,
    /// min_adc_integral, max_adc_integral, max_channel_span, max_time_span.
    /// An absent or null key keeps the default.
    TCMakerBurst::Config burst_config_from_json(const Configuration& cfg);

}  // namespace WireCell::Trig

#endif
