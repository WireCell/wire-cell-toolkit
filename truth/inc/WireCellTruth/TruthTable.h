/** TruthTable: the per-event particle tables of the truth file.

    The truth file truth_<ident>.npz is written by
    truth/python/dump_truth_table.py inside Stage A (LArSoft) and
    holds, per event, the tables of Truth::Schema:

      truth_nu   (Nn, 12) double   one row per generator interaction
      truth_pf   (Np, 25) double   one row per G4 particle, all kept
      sed_ids    (Ns, 2)  int32    per SimEnergyDeposit: trackid, orig_trackid
      sed_trackq (Nt, 2)  double   per trackid: visible electrons (IonAndScint)

    plus the string tables columns_nu / columns_pf / process_names /
    gen_labels as (K, L) uint8 ASCII matrices (cnpy reads no unicode).

    This class loads one file (cnpy: stored or deflated zip, not
    zip64), checks the column names against Schema and serves the
    rows; as_tensors() gives the two 2-D double tensors with
    md.datatype truth_nu / truth_pf and md.columns, as the SBND
    attacher wrote them, so that SbndPrMagnifyTrackingVisitor reads
    them unchanged.

    wcp-porting-img pr_ml/docs/19 sec 3.4.
 */

#ifndef WIRECELLTRUTH_TRUTHTABLE
#define WIRECELLTRUTH_TRUTHTABLE

#include "WireCellIface/ITensor.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Logging.h"

#include <map>
#include <string>
#include <vector>

namespace WireCell::Truth {

    struct TruthTable {
        /// Row-major tables; columns per Schema.
        std::vector<std::vector<double>> nu;
        std::vector<std::vector<double>> pf;
        /// Per SED row: (trackid, orig_trackid).  Empty when the file has none.
        std::vector<std::pair<int,int>> sed_ids;
        /// Per trackid: visible electrons from Stage A.  Empty when absent.
        std::map<int, double> sed_trackq;
        /// Generator labels, in nu_row order of the non-neutrino rows.
        std::vector<std::string> gen_labels;
        /// Where the table came from: "file", "depos" or "".
        std::string source;
        /// Free-form provenance from the file's "meta" member (JSON text), if any.
        Configuration meta;

        bool empty() const { return nu.empty() && pf.empty(); }

        /// Load from a truth file.  Throws IOError / ValueError on a
        /// missing file or a column-name mismatch with Schema.
        static TruthTable load(const std::string& fname, Log::logptr_t log = nullptr);

        /// Fallback with no particle table: one pf row per distinct
        /// depo key (prior()->id() when id_from_prior) with pdg,
        /// first time and position, no parentage, process unknown;
        /// one nu row (pdg 0, nu_idx -1).  source = "depos".
        static TruthTable from_depos(const IDepo::vector& depos, bool id_from_prior = true);

        /// pf row index of a trackid, or -1.
        int pf_row(int trackid) const;
        /// Column value helpers (throw KeyError on a bad name).
        double pf_value(size_t row, const std::string& col) const;
        double nu_value(size_t row, const std::string& col) const;

        /// (trackid, orig_trackid) of a depo id in Stage A index mode;
        /// (id, id) when sed_ids is empty or the id is out of range.
        std::pair<int,int> sed_id(int depo_id) const;

        /// The two tensors.  Datapaths "<prefix>/nu" and "<prefix>/pf";
        /// md: datatype, datapath, ident, columns.
        ITensor::vector as_tensors(int ident, const std::string& prefix) const;

      private:
        mutable std::map<int,int> m_pf_rows; // trackid -> row, built lazily
    };

}

#endif
