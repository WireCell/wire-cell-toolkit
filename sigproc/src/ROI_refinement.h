#ifndef WIRECELLSIGPROC_ROIREFINEMENT
#define WIRECELLSIGPROC_ROIREFINEMENT

#include "SignalROI.h"
#include "ROI_formation.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/RayGrid.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Logging.h"
#include "WireCellIface/IAnodePlane.h"

#include <vector>
#include <map>

namespace WireCell {
    namespace SigProc {

        class ROI_refinement {
           public:
            ROI_refinement(Waveform::ChannelMaskMap& cmm, int nwire_u, int nwire_v, int nwire_w, float th_factor = 3.0,
                           float fake_signal_low_th = 500, float fake_signal_high_th = 1000,
                           float fake_signal_low_th_ind_factor = 1.0, float fake_signal_high_th_ind_factor = 1.0,
                           int pad = 5, int break_roi_loop = 2, float th_peak = 3.0, float sep_peak = 6.0,
                           float low_peak_sep_threshold_pre = 1200, int max_npeaks = 200, float sigma = 2,
                           float th_percent = 0.1, bool isWrapped = false);
            ~ROI_refinement();

            void Clear();

            // initialize the ROIs
            void load_data(int plane, const Array::array_xxf& r_data, ROI_formation& roi_form);
            void refine_data(int plane, ROI_formation& roi_form);
            void refine_data_debug_mode(int plane, ROI_formation& roi_form, const std::string& cmd);

            void MP3ROI(const int plane, const IAnodePlane::pointer anode, const IAnodeFace::pointer face,
                        const std::map<int, int>& map_ch, ROI_formation& roi_form, const double mp_th1 = 0.,
                        const double mp_th2 = 0., const int tick_resolution = 10, const int wire_resolution = 2,
                        const int nbounds_layers = 2, const std::vector<int> iplane2layer={0,1,2});

            void MP2ROI(const int plane, const IAnodePlane::pointer anode, const IAnodeFace::pointer face,
                        const std::map<int, int>& map_ch, ROI_formation& roi_form, const double mp_th1 = 0.,
                        const double mp_th2 = 0., const int tick_resolution = 10, const int wire_resolution = 2,
                        const int nbounds_layers = 2, const std::vector<int> iplane2layer={0,1,2});

            typedef std::multimap<std::pair<int, int>, std::pair<int, int>> MapMPROI;
            const MapMPROI& get_mp2_rois() const { return mp_rois; }
            MapMPROI get_mp2_rois() { return mp_rois; }
            const MapMPROI& get_mp3_rois() const { return proteced_rois; }
            MapMPROI& get_mp3_rois() { return proteced_rois; }

            void CleanUpROIs(int plane);
            void generate_merge_ROIs(int plane);
            void CheckROIs(int plane, ROI_formation& roi_form);

            void CleanUpCollectionROIs();
            void CleanUpInductionROIs(int plane);
            void ShrinkROIs(int plane, ROI_formation& roi_form);
            void ShrinkROI(SignalROI* roi, ROI_formation& roi_form);

            void BreakROIs(int plane, ROI_formation& roi_form);
            void BreakROI(SignalROI* roi, float rms);
            void BreakROI1(SignalROI* roi);

            void ExtendROIs(int plane);

            void TestROIs();

            void apply_roi(int plane, Array::array_xxf& r_data);

            SignalROIChList& get_u_rois() { return rois_u_loose; };
            SignalROIChList& get_v_rois() { return rois_v_loose; };
            SignalROIChList& get_w_rois() { return rois_w_tight; };
            SignalROIChList& get_rois_by_plane(int plane);

           private:
            int nwire_u;
            int nwire_v;
            int nwire_w;

            float th_factor;
            float fake_signal_low_th;
            float fake_signal_high_th;
            float fake_signal_low_th_ind_factor;
            float fake_signal_high_th_ind_factor;
            int pad;
            int break_roi_loop;
            float th_peak;
            float sep_peak;
            float low_peak_sep_threshold_pre;

            // peak finding ones
            int max_npeaks;
            float sigma;
            float th_percent;

            void unlink(SignalROI* prev_roi, SignalROI* next_roi);
            void link(SignalROI* prev_roi, SignalROI* next_roi);

            std::map<int, std::vector<std::pair<int, int>>> bad_ch_map;

            SignalROIChList rois_u_tight;
            SignalROIChList rois_v_tight;
            SignalROIChList rois_w_tight;

            SignalROIChList rois_u_loose;
            SignalROIChList rois_v_loose;

            MapMPROI proteced_rois;  // using chid and start_bin as id
            MapMPROI mp_rois;        // using chid and start_bin as id

            SignalROIMap front_rois;
            SignalROIMap back_rois;
            SignalROIMap contained_rois;

            Log::logptr_t log;

            bool isWrapped;
        };
    }  // namespace SigProc
}  // namespace WireCell
#endif
// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
