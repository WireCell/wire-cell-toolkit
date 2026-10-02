/** FdvdLowEQLMatching: FD-VD low-energy (solar) charge-light matching.
 *
 * The frozen python matcher of fdvd_sim docs 06-10 and 15 (purity-first,
 * strict one-to-one; numeric core in FdvdLowE.h) as a WCT fan-in node.
 *
 * Inputs (ITensorSet ports):
 *   0  clustering point-cloud tree (MultiAlgBlobClustering output; read at
 *      inpath/live, inpath "%d" formatted with the set ident)
 *   1  blob table (FdvdBlobTable): light-prediction blob centres and the join
 *      key to the tree's blob scalars
 *   2  flashes: the OpFlashFinder tensor set ("opflash" [nflash, 1 + 184],
 *      column 0 time [ns])
 *   3  optional, multiplicity 4: drift values from elsewhere ("drift" [n, 3] =
 *      representative blob order, mu, sigma [cm]), e.g. a file
 * Drift: with a "drift" configuration block the matcher runs the drift
 * regressor itself (FdvdDriftRegressor, multiplicity 3); with neither, no
 * drift veto can cut (mu NaN, as python for clusters without a crop).
 *
 * Clusters: the tree's clusters with blob charge sum Q >= qmin_cluster whose
 * shift table exists (ql_m2m_proto.collect_one + ql08.cl_one without the
 * truth-only MARLEY-bundle inclusion: such a cluster is below Qc = 100 ke
 * and cannot take part in a decision, fdvd_sim doc 16 sec 2).
 *
 * Output: one ITensorSet (ident of port 0) with tensors
 *   clusters   f8 [ncl, 10]  tree index, rep blob, nblob, Q, x_app, t_lo, t_hi, dt, mu, sigma
 *   membership f8 [nb, 2]    tree cluster index, blob-table order  (every tree blob)
 *   groups     f8 [ng, 5 + 184]  window index, group index, t [us], tot PE, lit OpDets, PE[184]
 *   pairs      f8 [np, 6 + nspec]  k, g, r, ks, dc, npp, drift-pass per spec
 *   decisions  f8 [nd, 4]    arm, spec, k, g   (k: clusters row, g: groups window index)
 *   drift      f8 [n, 3]    rep, mu, sigma of every regressed cluster (drift block only)
 *   crops      f4 [n, 256, 1024]  the regressor's own crops (drift.dump_crops only)
 *   tables     f4 [sum nstep, 184] + table_index f8 [ncl, 2] + the input blobs table  (dump_tables only)
 * The set metadata names the arms and specs and counts join mismatches.
 */
#ifndef WIRECELL_MATCH_FDVDLOWEQLMATCHING
#define WIRECELL_MATCH_FDVDLOWEQLMATCHING

#include "WireCellIface/ITensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellMatch/FdvdLowE.h"
#include "WireCellMatch/FdvdDriftRegressor.h"

#include <memory>

namespace WireCell::Match {

    class FdvdLowEQLMatching : public Aux::Logger, public ITensorSetFanin, public IConfigurable {
      public:
        FdvdLowEQLMatching();
        virtual ~FdvdLowEQLMatching();

        virtual std::vector<std::string> input_types();
        virtual bool operator()(const input_vector& invec, output_pointer& out);

        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

      private:
        size_t m_multiplicity{3};
        std::string m_inpath{"pointtrees/%d"};
        std::string m_library{""};       // PhotonLibraryModel meta JSON (fdvd-photlib-vis-comb-10cm.json)
        std::string m_calibration{""};   // export_lowe_cal.py JSON
        std::string m_geom_file{""};     // OpDet geometry JSON {"opdets": [{opdet, x, y, z} mm]}
        std::vector<double> m_veto_k{0.0, 3.0};
        std::vector<std::string> m_spec_names{"none", "veto3"};
        bool m_dump_tables{false};
        FdvdLowE::Constants m_C;
        FdvdLowE::Calibration m_cal;
        std::vector<FdvdLowE::Arm> m_arms;
        std::vector<std::array<double, 3>> m_pd_pos;
        std::unique_ptr<PhotonLibraryModel> m_lib;
        std::unique_ptr<FdvdDriftRegressor> m_drift;
        size_t m_count{0};
    };

}  // namespace WireCell::Match

#endif
