/** FdvdLowERootWriter: the FD-VD low-energy reconstruction result as a ROOT file.
 *
 * A pass-through filter on the tensor set of Match::FdvdLowEQLMatching (fdvd_sim doc 16): the set is
 * forwarded unchanged (so a TensorFileSink can still archive it) and written as flat TTrees, one
 * set of rows per event (branch "event" = the set ident):
 *   T_meta     one entry per job: arm and spec names (comma separated), producer
 *   T_event    counts: tree clusters / blobs, unjoined blobs, flashes, stored and in-window flash groups,
 *              Q-L clusters, pairs, regressed clusters
 *   T_cluster  per Q-L cluster (Q >= 50 ke with a light-prediction table): tree index, representative blob,
 *              nblob, Q, x_app, drift window [t_lo, t_hi], drift regressor mu / sigma (NaN if none), and per
 *              arm x spec (vectors, index arm * nspec + spec) the matched flash group (-1), its time t0 [us]
 *              and the drift-corrected x = x_app + vcm t0 [cm]
 *   T_flash    per in-window flash group: window index, group index, time [us], total PE, lit OpDets, PE[184]
 *   T_pair     per (cluster, flash group) pair in the drift window: k, g, r, ks, dc, npp, pass[spec];
 *              which pairs: config "pairs" = "qc" (default: cluster Q >= pair_qmin and group >= pair_npd lit
 *              OpDets, the pairs that can enter a decision, as the python t1 record keeps), "all" (~130 k
 *              per event) or "none"
 *   T_match    per decision: arm, spec, k, g, tree index, representative blob, t0 [us], drift-corrected x [cm]
 *   T_drift    per regressed cluster: representative blob, mu, sigma [cm]
 * k indexes T_cluster rows of the same event, g T_flash rows (window index).  The file is written at EOS.
 *
 * Configuration: output_filename, vcm (cm/us, default the FD-VD 0.160563), pairs, pair_qmin (100e3),
 * pair_npd (3), event (>= 0: the "event"
 * value written instead of the set ident, e.g. for one-event jobs whose inputs carry a fixed ident).
 */
#ifndef WIRECELLROOT_FDVDLOWEROOTWRITER
#define WIRECELLROOT_FDVDLOWEROOTWRITER

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <memory>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace WireCell::Root {

    class FdvdLowERootWriter : public Aux::Logger, public ITensorSetFilter, public IConfigurable {
      public:
        FdvdLowERootWriter();
        virtual ~FdvdLowERootWriter();

        virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

      private:
        void open(const ITensorSet::pointer& in);
        void close();

        std::string m_filename{"fdvd-lowe.root"};
        double m_vcm{0.160563};
        int m_event{-1};
        std::string m_pairs{"qc"};
        double m_pair_qmin{100e3};
        int m_pair_npd{3};
        TFile* m_file{nullptr};
        TTree *m_tmeta{nullptr}, *m_tevent{nullptr}, *m_tcluster{nullptr}, *m_tflash{nullptr}, *m_tpair{nullptr},
            *m_tmatch{nullptr}, *m_tdrift{nullptr};
        size_t m_narm{0}, m_nspec{0};
        struct Buf;                      // branch buffers (FdvdLowERootWriter.cxx)
        std::unique_ptr<Buf> m_buf;
        size_t m_count{0};
    };

}  // namespace WireCell::Root

#endif
