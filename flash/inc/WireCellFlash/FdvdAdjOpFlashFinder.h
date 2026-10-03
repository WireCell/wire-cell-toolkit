/** AdjOpHits flash finder for FD-VD (DUNE far detector vertical drift).
 *
 * A port of duneopdet LowEPDSUtils AdjOpHitsUtils::CalcAdjOpHits (v10_26_00d00,
 * AdjOpHitsUtils.cc:262-494; the plane assignment GetOpHitPlane / CheckPlane
 * :559-613; the flash vector essentials of MakeFlashVector :29-216), by way of
 * its python port fdvd_sim/stage19/adjophits.py (wcp
 * WireCell/wcp-porting-validation, checked against the real SolarOpFlash
 * module on 50/50 events, fdvd_sim doc 19) and the nanosecond form
 * fdvd_sim/stage21/adj21.py, which this component reproduces bit for bit
 * (fdvd_sim doc 21 gate A1).
 *
 * The C++ quirks of CalcAdjOpHits are kept on purpose:
 *  - the backward scan never reaches the first hit of the time-sorted list
 *    (`it4 != begin`);
 *  - a brighter candidate demotes the seed and un-claims every hit collected
 *    so far, including hits an earlier flash had claimed;
 *  - two hits on no plane (-1) count as the same plane.
 *
 * Input: ITensorSet with an "ophits" tensor (f8 [nhit, 9], the FdvdOpHitFinder
 * / OpFlashFinder column order, times in WCT ns): 0 OpChannel, 1 peak_time,
 * 5 PE, 6 start_time.  OpChannel -> OpDet through channel_map_file.
 *
 * Output: ITensorSet holding
 *   "opflash"        f8 [nflash, 1 + nchan]: time [ns] = time of the max-PE
 *                    member (time_var), then the member PE summed per OpDet;
 *                    the format FdvdLowEQLMatching and OpFlashFinder use
 *   "flash_summary"  f8 [nflash, 8]: index, total PE, hot-vertex y, z [cm],
 *                    plane, members, PE-weighted time [ns], trigger hit row
 *   "adj_membership" f8 [nmember, 2]: (flash, ophits row); with
 *                    hit_duplicates a hit may sit in several flashes
 *   "ophits"         the input hits, passed through.
 */
#ifndef WIRECELLFLASH_FDVDADJOPFLASHFINDER
#define WIRECELLFLASH_FDVDADJOPFLASHFINDER

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <array>
#include <map>
#include <string>
#include <vector>

namespace WireCell {
    namespace Flash {

        // CalcAdjOpHits parameters.  Defaults = the official FD-VD SolarOpFlash
        // setting solar_opflash_dunevd10kt_1x8x14_3view_30deg (SolarOpFlash.fcl:30-41).
        struct FdvdAdjParams {
            bool start_time{true};          // OpFlashAlgoTime "StartTime" (else PeakTime)
            double min_time_us{0.030};      // backward window
            double max_time_us{0.030};      // forward window
            double radius_cm{500.0};        // OpDet centre distance
            int nhit{3};                    // minimum members
            double pe{1.5};                 // member PE threshold
            double trigger_pe{1.5};         // seed PE threshold
            double hot_threshold{0.3};      // hot vertex: members >= hot_threshold * max PE
            bool hit_duplicates{true};      // a claimed hit may join another flash
            // X-ARAPUCA plane positions [cm] (GetOpHitPlane) and the CheckPlane buffer
            double xa_cathode_x{-327.5};
            double xa_membrane_y{743.302};
            double xa_final_cap_z{2188.38};
            double xa_start_cap_z{-96.5};
            double plane_buffer{0.1};
        };

        // One flash: members as row indices into the time-sorted input order
        // of the caller's arrays (element 0 = the trigger hit).
        using FdvdAdjCluster = std::vector<size_t>;

        // CalcAdjOpHits on hit arrays (time [ns], PE, OpDet); dist[a][b] = OpDet
        // centre distance [cm], plane[od] = GetOpHitPlane.
        std::vector<FdvdAdjCluster> fdvd_calc_adj_ophits(const std::vector<double>& t_ns, const std::vector<double>& pe,
                                                         const std::vector<int>& opdet,
                                                         const std::vector<std::vector<double>>& dist,
                                                         const std::vector<int>& plane, const FdvdAdjParams& par);

        // GetOpHitPlane per OpDet centre [cm]: 0 cathode, 1 membrane +y, 2 membrane -y,
        // 3 final cap, 4 start cap, -1 none.
        int fdvd_adj_plane(double x, double y, double z, const FdvdAdjParams& par);

        class FdvdAdjOpFlashFinder : public Aux::Logger, public ITensorSetFilter, public IConfigurable {
          public:
            FdvdAdjOpFlashFinder();
            virtual ~FdvdAdjOpFlashFinder();

            virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            FdvdAdjParams m_par;
            int m_nchan{184};
            std::string m_geom_file{""};
            std::string m_channel_map_file{""};
            std::map<int, int> m_chmap;                  // OpChannel -> OpDet
            std::vector<std::array<double, 3>> m_pos;    // OpDet centres [cm]
            std::vector<std::vector<double>> m_dist;     // centre distances [cm]
            std::vector<int> m_plane;                    // plane per OpDet
            size_t m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
