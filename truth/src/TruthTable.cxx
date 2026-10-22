#include "WireCellTruth/TruthTable.h"
#include "WireCellTruth/Schema.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/cnpy.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"

#include <algorithm>
#include <map>

using namespace WireCell;

namespace {

    // A (K, L) uint8 ASCII matrix, zero padded, as strings.
    std::vector<std::string> load_strings(const cnpy::npz_t& npz, const std::string& name)
    {
        std::vector<std::string> out;
        auto it = npz.find(name);
        if (it == npz.end()) return out;
        const auto& arr = it->second;
        if (arr.shape.size() != 2 || arr.word_size != 1) {
            THROW(ValueError() << errmsg{"truth file member " + name + " is not a (K,L) uint8 matrix"});
        }
        const size_t nrow = arr.shape[0], ncol = arr.shape[1];
        const char* d = arr.data<char>();
        for (size_t r = 0; r < nrow; ++r) {
            std::string s(d + r * ncol, ncol);
            s.erase(std::find(s.begin(), s.end(), '\0'), s.end());
            out.push_back(s);
        }
        return out;
    }

    template<typename T>
    std::vector<std::vector<double>> load_table(const cnpy::npz_t& npz, const std::string& name, size_t ncols)
    {
        std::vector<std::vector<double>> out;
        auto it = npz.find(name);
        if (it == npz.end()) return out;
        const auto& arr = it->second;
        if (arr.shape.size() != 2) {
            THROW(ValueError() << errmsg{"truth file member " + name + " is not 2-D"});
        }
        if (arr.word_size != sizeof(T)) {
            THROW(ValueError() << errmsg{"truth file member " + name + " has an unexpected word size"});
        }
        const size_t nrow = arr.shape[0], ncol = arr.shape[1];
        if (ncols && ncol != ncols) {
            THROW(ValueError() << errmsg{"truth file member " + name + " has " + std::to_string(ncol) +
                                         " columns, schema has " + std::to_string(ncols)});
        }
        const T* d = arr.data<T>();
        out.resize(nrow, std::vector<double>(ncol, 0.0));
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                // numpy writes C order; cnpy reports fortran_order when not
                out[r][c] = arr.fortran_order ? (double)d[c * nrow + r] : (double)d[r * ncol + c];
            }
        }
        return out;
    }

    void check_columns(const std::vector<std::string>& got, const std::vector<std::string>& want,
                       const std::string& what)
    {
        if (got.empty()) return;  // an old file without the name table: trust the shape
        if (got != want) {
            std::string g, w;
            for (const auto& s : got) g += s + ",";
            for (const auto& s : want) w += s + ",";
            THROW(ValueError() << errmsg{"truth file " + what + " columns [" + g + "] differ from schema [" + w + "]"});
        }
    }

    Configuration columns_md(const std::vector<std::string>& cols)
    {
        Configuration jcols = Json::arrayValue;
        for (const auto& c : cols) jcols.append(c);
        return jcols;
    }
}

Truth::TruthTable Truth::TruthTable::load(const std::string& fname, Log::logptr_t log)
{
    TruthTable tt;
    cnpy::npz_t npz;
    try {
        npz = cnpy::npz_load(fname);
    }
    catch (std::runtime_error& err) {
        THROW(IOError() << errmsg{"failed to load truth file " + fname + ": " + err.what()});
    }
    check_columns(load_strings(npz, "columns_nu"), truth_nu_columns(), "truth_nu");
    check_columns(load_strings(npz, "columns_pf"), truth_pf_columns(), "truth_pf");
    tt.nu = load_table<double>(npz, truth_file_nu, truth_nu_columns().size());
    tt.pf = load_table<double>(npz, truth_file_pf, truth_pf_columns().size());

    const auto pnames = load_strings(npz, "process_names");
    if (!pnames.empty() && pnames != process_names()) {
        // The file's codes were assigned from its own list; remap to ours.
        const int pcol = column_index(truth_pf_columns(), "process");
        size_t nremap = 0;
        for (auto& row : tt.pf) {
            const int code = (int)row[pcol];
            if (code < 0 || code >= (int)pnames.size()) continue;
            const int ours = process_code(pnames[code]);
            if (ours != code) { row[pcol] = ours; ++nremap; }
        }
        if (log) log->warn("truth file {} process list differs from schema: {} codes remapped", fname, nremap);
    }

    {
        auto it = npz.find(truth_file_sed_ids);
        if (it != npz.end()) {
            const auto& arr = it->second;
            if (arr.shape.size() != 2 || arr.shape[1] != 2 || arr.word_size != 4) {
                THROW(ValueError() << errmsg{"truth file sed_ids must be (N,2) int32"});
            }
            const int32_t* d = arr.data<int32_t>();
            const size_t n = arr.shape[0];
            tt.sed_ids.reserve(n);
            for (size_t r = 0; r < n; ++r) {
                if (arr.fortran_order) tt.sed_ids.emplace_back(d[r], d[n + r]);
                else tt.sed_ids.emplace_back(d[2 * r], d[2 * r + 1]);
            }
        }
    }
    {
        auto it = npz.find("sed_energy");
        if (it != npz.end() && it->second.word_size == 4 && it->second.shape.size() == 1) {
            const float* d = it->second.data<float>();
            tt.sed_energy.assign(d, d + it->second.shape[0]);
        }
    }
    for (const auto& row : load_table<double>(npz, truth_file_sed_trackq, 2)) {
        tt.sed_trackq[(int)row[0]] = row[1];
    }
    tt.gen_labels = load_strings(npz, "gen_labels");
    {
        const auto meta = load_strings(npz, "meta");
        if (!meta.empty()) {
            try { tt.meta = Persist::loads(meta[0]); }
            catch (...) { if (log) log->warn("truth file {}: meta is not JSON", fname); }
        }
    }
    tt.source = "file";
    if (log) {
        log->debug("loaded truth file {}: {} nu rows, {} pf rows, {} sed ids, {} track charges",
                   fname, tt.nu.size(), tt.pf.size(), tt.sed_ids.size(), tt.sed_trackq.size());
    }
    return tt;
}

Truth::TruthTable Truth::TruthTable::from_depos(const IDepo::vector& depos, bool id_from_prior)
{
    TruthTable tt;
    struct Seen { int pdg; double t; Point pos; double q; };
    std::map<int, Seen> first;
    for (const auto& d : depos) {
        if (!d) continue;
        IDepo::pointer key = d;
        if (id_from_prior) {
            while (key->prior()) key = key->prior();
        }
        const int id = key->id();
        auto it = first.find(id);
        if (it == first.end()) {
            first[id] = Seen{key->pdg(), key->time(), key->pos(), std::abs(d->charge())};
        }
        else {
            if (key->time() < it->second.t) {
                it->second.t = key->time();
                it->second.pos = key->pos();
            }
            it->second.q += std::abs(d->charge());
        }
    }
    const auto& nuc = truth_nu_columns();
    std::vector<double> nurow(nuc.size(), 0.0);
    nurow[column_index(nuc, "nu_idx")] = -1;
    tt.nu.push_back(nurow);

    const auto& pfc = truth_pf_columns();
    for (const auto& [id, s] : first) {
        std::vector<double> row(pfc.size(), 0.0);
        row[column_index(pfc, "nu_row")] = 0;
        row[column_index(pfc, "trackid")] = id;
        row[column_index(pfc, "parent_trackid")] = 0;
        row[column_index(pfc, "mother_trackid")] = 0;
        row[column_index(pfc, "pdg")] = s.pdg;
        row[column_index(pfc, "process")] = process_unknown;
        row[column_index(pfc, "start_x")] = s.pos.x();
        row[column_index(pfc, "start_y")] = s.pos.y();
        row[column_index(pfc, "start_z")] = s.pos.z();
        row[column_index(pfc, "start_t")] = s.t;
        row[column_index(pfc, "end_x")] = s.pos.x();
        row[column_index(pfc, "end_y")] = s.pos.y();
        row[column_index(pfc, "end_z")] = s.pos.z();
        row[column_index(pfc, "end_t")] = s.t;
        row[column_index(pfc, "orig_trackid")] = id;
        row[column_index(pfc, "gen_label")] = -1;
        row[column_index(pfc, "pf_node")] = -1;
        row[column_index(pfc, "cls")] = (double)(int)Cls::unknown;
        row[column_index(pfc, "merged_into")] = -1;
        row[column_index(pfc, "nvis_q")] = s.q;
        tt.pf.push_back(row);
    }
    tt.source = "depos";
    return tt;
}

int Truth::TruthTable::pf_row(int trackid) const
{
    if (m_pf_rows.empty() && !pf.empty()) {
        const int tcol = column_index(truth_pf_columns(), "trackid");
        for (size_t r = 0; r < pf.size(); ++r) {
            m_pf_rows.emplace((int)pf[r][tcol], (int)r);  // first wins
        }
    }
    auto it = m_pf_rows.find(trackid);
    return it == m_pf_rows.end() ? -1 : it->second;
}

double Truth::TruthTable::pf_value(size_t row, const std::string& col) const
{
    const int c = column_index(truth_pf_columns(), col);
    if (c < 0) THROW(KeyError() << errmsg{"no truth_pf column " + col});
    return pf.at(row).at(c);
}

double Truth::TruthTable::nu_value(size_t row, const std::string& col) const
{
    const int c = column_index(truth_nu_columns(), col);
    if (c < 0) THROW(KeyError() << errmsg{"no truth_nu column " + col});
    return nu.at(row).at(c);
}

std::pair<int,int> Truth::TruthTable::sed_id(int depo_id) const
{
    if (depo_id < 0 || (size_t)depo_id >= sed_ids.size()) return {depo_id, depo_id};
    return sed_ids[depo_id];
}

ITensor::vector Truth::TruthTable::as_tensors(int ident, const std::string& prefix) const
{
    ITensor::vector out;
    auto make = [&](const std::vector<std::vector<double>>& table, const std::vector<std::string>& cols,
                    const std::string& datatype, const std::string& datapath) {
        std::vector<double> flat;
        flat.reserve(table.size() * cols.size());
        for (const auto& row : table) flat.insert(flat.end(), row.begin(), row.end());
        Configuration md;
        md["datatype"] = datatype;
        md["datapath"] = datapath;
        md["ident"] = ident;
        md["columns"] = columns_md(cols);
        md["schema_version"] = schema_version;
        ITensor::shape_t shape{table.size(), cols.size()};
        out.push_back(std::make_shared<Aux::SimpleTensor>(shape, flat.data(), md));
    };
    make(nu, truth_nu_columns(), truth_nu_datatype, prefix + "/nu");
    make(pf, truth_pf_columns(), truth_pf_datatype, prefix + "/pf");
    return out;
}
