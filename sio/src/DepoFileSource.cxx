#include "WireCellSio/DepoFileSource.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"

#include "WireCellAux/DepoTools.h"

#include <map>
#include <mutex>

WIRECELL_FACTORY(DepoFileSource, WireCell::Sio::DepoFileSource,
                 WireCell::INamed,
                 WireCell::IDepoSetSource, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Stream;
using namespace WireCell::String;

// ident, type, tag, ext
struct filemd_t {
    std::string fname{""};
    bool okay{false};
    int ident{0};
    std::string type{""}, tag{""}, ext{""};
    // "depo", "data|info|truth", "npy
};

// A member read ahead of its use.
struct WireCell::Sio::DepoFilePending {
    pigenc::File pig;
    filemd_t fmd;
};

Sio::DepoFileSource::DepoFileSource()
    : Aux::Logger("DepoFileSource", "io")
{
}

Sio::DepoFileSource::~DepoFileSource()
{
}

WireCell::Configuration Sio::DepoFileSource::default_configuration() const
{
    Configuration cfg;
    // Input tar stream
    cfg["inname"] = m_inname;

    // This will be multiplied to the "charge" of each depo.  Eg, it
    // can be usfully set to -1 to accomodate diffierent conventions
    // for electron charge sign.
    cfg["scale"] = m_scale;

    // pr_ml doc 33: one loaded copy for all sources of the same file
    // and scale in this process (see the header).  Default false.
    cfg["share_loaded"] = m_share;
    return cfg;
}

void Sio::DepoFileSource::configure(const WireCell::Configuration& cfg)
{
    m_inname = get(cfg, "inname", m_inname);
    m_scale = get(cfg, "scale", m_scale);
    m_share = get(cfg, "share_loaded", m_share);

    m_in.clear();
    m_pending.reset();
    input_filters(m_in, m_inname);
    if (m_in.empty()) {
        THROW(ValueError() << errmsg{"DepoFileSource: unsupported inname: " + m_inname});
    }

    log->debug("reading {} with scale={}", m_inname, m_scale);
}

static
filemd_t parse_filename(const std::string& fname)
{
    if (fname.empty()) { return {fname}; }
    auto parts = split(fname, "_");
    if (parts.size() != 3) { return {fname}; }
    auto rparts = split(parts[2], ".");
    int ident = std::atoi(rparts[0].c_str());
    return filemd_t{fname, true, ident, parts[0], parts[1], rparts[1]};
}

static
filemd_t read_pig(std::istream& si, pigenc::File& pig)
{
    std::string fname{""};
    size_t fsize{0};
    custard::read(si, fname, fsize);
    pig.read(si);
    return parse_filename(fname);
}

// Read one file set from tar stream and save to frame if tag matches
IDepoSet::pointer Sio::DepoFileSource::next()
{
    // Full eager load of next two files, the first possibly already
    // read ahead by the previous call while looking for a sidecar.
    pigenc::File pig1, pig2;
    filemd_t fmd1;
    if (m_pending) {
        pig1 = std::move(m_pending->pig);
        fmd1 = m_pending->fmd;
        m_pending.reset();
    }
    else {
        fmd1 = read_pig(m_in, pig1);
        if (m_in.eof()) {
            log->debug("call={}, read1 depo stream EOF with file={}", m_count, m_inname);
            return nullptr;
        }
    }
    auto fmd2 = read_pig(m_in, pig2);
    if (m_in.eof()) {
        log->debug("call={}, read2 depo stream EOF with file={}", m_count, m_inname);
        return nullptr;
    }
    if (!m_in) {
        log->error("call={}, depo read fail with file={}", m_count, m_inname);
        return nullptr;
    }
    if (!fmd1.okay) {
        log->error("call={}, depo read fail to parse npy file name {} in file={}", m_count, fmd1.fname, m_inname);
        return nullptr;
    }
    if (!fmd2.okay) {
        log->error("call={}, depo read fail to parse npy file name {} in file={}", m_count, fmd2.fname, m_inname);
        return nullptr;
    }
    if (fmd1.ident != fmd2.ident) {
        log->error("call={}, ident mismatch {} != {} in file={}",
                   m_count, m_count, fmd1.ident, fmd2.ident);
        return nullptr;
    }
    const int ident = fmd1.ident;

    // Allow free order of data/info w/in the pair.
    pigenc::File* dpig=&pig1;
    pigenc::File* ipig=&pig2;
    if (fmd2.tag == "data") {
        dpig = &pig2;
        ipig = &pig1;
    }
    /// Note, this emits have non-ascii chars from Numpy header
    // log->debug("call={}, data header: {}", m_count, dpig->header().str());
    // log->debug("call={}, info header: {}", m_count, ipig->header().str());

    // Optional third member: the truth sidecar.  Read one member
    // ahead; if it is not ours, keep it for the next call.
    using array_xxdrw = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    array_xxdrw tarr;
    bool have_truth = false;
    {
        auto pend = std::make_unique<DepoFilePending>();
        pend->fmd = read_pig(m_in, pend->pig);
        if (m_in.eof() || !m_in) {
            // last set of the stream, no sidecar; the next call will
            // see the same EOF.
        }
        else if (pend->fmd.okay && pend->fmd.tag == "truth" && pend->fmd.ident == ident) {
            if (pigenc::eigen::load(pend->pig, tarr)) {
                have_truth = true;
            }
            else {
                log->warn("call={}, failed to load truth sidecar {} in file={}, ignoring",
                          m_count, pend->fmd.fname, m_inname);
            }
        }
        else {
            m_pending = std::move(pend);
        }
    }

    using array_xxfrw = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    array_xxfrw darr;
    using array_xxirw = Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    array_xxirw iarr;

    bool ok = false;

    ok = pigenc::eigen::load(*dpig, darr);
    if (!ok) {
        log->error("call={}, short load failed data array file={}",
                   m_count, m_inname);
        return nullptr;
    }

    ok = pigenc::eigen::load(*ipig, iarr);
    if (!ok) {
        log->error("call={}, short load failed for info array file={}",
                   m_count, m_inname);
        return nullptr;
    }

    if (darr.rows() != iarr.rows() and darr.cols() == iarr.cols()) {
        array_xxfrw dtmp = darr.transpose();
        darr = dtmp;
        array_xxirw itmp = iarr.transpose();
        iarr = itmp;
        log->warn("call={}, array transpose detected, correcting to: data=({},{}) info=({},{})",
                  m_count, darr.rows(), darr.cols(), iarr.rows(), iarr.cols());
    }

    if (darr.rows() != iarr.rows()) {
        log->error("call={}, array row mismatch. transposed? data=({},{}) info=({},{})",
                   m_count, darr.rows(), darr.cols(), iarr.rows(), iarr.cols());
        return nullptr;
    }
    const size_t ndepos = darr.rows();

    const size_t ndata = 7;
    if (darr.cols() != ndata) {
        log->error("call={}, data column mismatch got={} want={}",
                   m_count, darr.cols(), ndata);
        return nullptr;
    }

    const size_t ninfo = 4;
    if (iarr.cols() != ninfo) {
        log->error("call={}, info column mismatch got={} want={}",
                   m_count, iarr.cols(), ninfo);
        return nullptr;
    }

    const size_t ntruth = 3;
    if (have_truth) {
        if ((size_t)tarr.rows() != ndepos || (size_t)tarr.cols() != ntruth) {
            log->warn("call={}, truth sidecar shape mismatch got=({},{}) want=({},{}), ignoring",
                      m_count, tarr.rows(), tarr.cols(), ndepos, ntruth);
            have_truth = false;
        }
    }

    // Files written before 2026-10-21 carry child == own row index on
    // every gen>0 row (see Aux::fill).  Rows are pushed depth-first,
    // so in that case the depo a prior belongs to is the previous row.
    size_t nprior_rows = 0, nself = 0;
    for (size_t ind=0; ind < ndepos; ++ind) {
        if (iarr(ind, 2) > 0) {
            ++nprior_rows;
            if ((size_t)iarr(ind, 3) == ind) {
                ++nself;
            }
        }
    }
    const bool legacy_child = nprior_rows > 0 && nself == nprior_rows;
    if (legacy_child) {
        log->warn("call={}, legacy self-referencing child indices on {} prior rows in file={}, using row order",
                  m_count, nprior_rows, m_inname);
    }

    // build list of depos in order given
    std::vector<std::shared_ptr<Aux::SimpleDepo>> all_sdepos;

    // First pass make all depos
    for (size_t ind=0; ind < ndepos; ++ind) {

        const double energy = have_truth ? tarr(ind, 0) : 1.0;
        auto sdepo = std::make_shared<Aux::SimpleDepo>(
            darr(ind, 0),        // t
            Point(darr(ind, 2),  // x
                  darr(ind, 3),  // y
                  darr(ind, 4)), // z
            darr(ind, 1)*m_scale,// q
            nullptr,             // prior
            darr(ind, 5),        // extent_long
            darr(ind, 6),        // extent_tran
            iarr(ind, 0),        // id
            iarr(ind, 1),        // pdg
            energy);

        all_sdepos.push_back(sdepo);
    }

    // second pass, save out the active and resolve the prior depos
    std::vector<std::shared_ptr<Aux::SimpleDepo>> sdepos;
    size_t npriors = 0;
    size_t npriors_missing = 0;
    for (size_t ind=0; ind < ndepos; ++ind) {

        const auto gen = iarr(ind, 2);
        if (!gen) {  // active
            sdepos.push_back(all_sdepos[ind]);
            continue;
        }

        // Prior: "other" is the row of the depo whose prior this is.
        size_t other = iarr(ind, 3);
        if (legacy_child) {
            other = ind - 1;    // ind > 0 as a gen>0 row follows its child
            if (iarr(other, 2) != gen - 1) {
                ++npriors_missing;
                continue;
            }
        }
        if (other >= all_sdepos.size() || other == ind) {
            ++npriors_missing;
        }
        else {
            ++npriors;
            all_sdepos[other]->set_prior(all_sdepos[ind]);
        }
    }
    if (npriors_missing) {
        log->warn("call={}, missing {} prior depos, active={} total={}",
                  m_count, npriors_missing, sdepos.size(), all_sdepos.size());
    }
    all_sdepos.clear();

    // convert to IDepos
    WireCell::IDepo::vector idepos;
    for (auto sdepo: sdepos) {
        idepos.push_back(sdepo);
    }
    sdepos.clear();

    log->debug("call={} loaded {} active, {} total depos from ident {} stream {} with {} priors, truth sidecar={}",
               m_count, idepos.size(), ndepos, ident, m_inname, npriors, have_truth);

    return std::make_shared<Aux::SimpleDepoSet>(ident, idepos);
}

namespace {
    // share_loaded: (inname, scale) -> the loaded sets, held only by the
    // sources still emitting them and by their consumers
    using shared_sets_t = std::vector<IDepoSet::pointer>;
    std::mutex g_share_mutex;
    std::map<std::pair<std::string, double>, std::weak_ptr<const shared_sets_t>> g_share;
}

bool Sio::DepoFileSource::operator()(IDepoSet::pointer& ds)
{
    ds = nullptr;
    if (m_share) {
        if (m_eos_sent) {
            ++m_count;
            return false;
        }
        if (!m_all) {
            std::lock_guard<std::mutex> lock(g_share_mutex);
            auto& slot = g_share[std::make_pair(m_inname, m_scale)];
            m_all = slot.lock();
            if (m_all) {
                log->debug("call={} sharing {} loaded set(s) of {}", m_count, m_all->size(), m_inname);
            }
            else {
                auto sets = std::make_shared<shared_sets_t>();
                while (auto one = next()) {
                    sets->push_back(one);
                }
                m_all = sets;
                slot = m_all;
                log->debug("call={} loaded {} set(s) of {} to share", m_count, m_all->size(), m_inname);
            }
        }
        if (m_index < m_all->size()) {
            ds = (*m_all)[m_index++];
        }
        else {
            m_all.reset();
            m_eos_sent = true;
            log->debug("EOS at call={}", m_count);
        }
        ++m_count;
        return true;
    }
    if (m_eos_sent) {
        ++m_count;
        return false;
    }
    ds = next();
    if (!ds) {
        m_eos_sent = true;
        log->debug("EOS at call={}", m_count);
    }
    ++m_count;
    return true;
}
