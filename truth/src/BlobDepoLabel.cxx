#include "WireCellTruth/BlobDepoLabel.h"
#include "WireCellTruth/TruthTable.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(BlobDepoLabel, WireCell::Truth::BlobDepoLabel,
                 WireCell::INamed,
                 WireCell::IBlobDepoFill, WireCell::IConfigurable)

using namespace WireCell;

Configuration Truth::assoc_default_configuration(const DepoBlobAssocParams& p)
{
    Configuration cfg;
    // As BlobDepoFill: added to depo times before associating with slices.
    cfg["time_offset"] = p.time_offset;
    // As BlobDepoFill: number of sigma of the depo Gaussian to consider.
    cfg["nsigma"] = p.nsigma;
    // As BlobDepoFill: drift speed converting extent_long to time.  Required.
    cfg["speed"] = p.speed;
    // As BlobDepoFill: the primary plane the depos are diced onto.
    cfg["pindex"] = p.pindex;
    // Track key and energy from the depo's prior() (the undrifted G4
    // depo) when it exists, else from the depo itself.
    cfg["id_from_prior"] = p.id_from_prior;
    // Extra widths added in quadrature to the depo extents, in time
    // (longitudinal) and distance (transverse), e.g. the signal
    // processing filter widths.  0 = off: the integration is then
    // exactly BlobDepoFill's and the filled cluster is bit-identical.
    cfg["extra_smear_long"] = p.extra_smear_long;
    cfg["extra_smear_tran"] = p.extra_smear_tran;
    // Per-plane transverse widths (collection across the primary wires,
    // induction across the other planes' wires; see DepoBlobAssocParams).
    // Written only when on, as point_zero_width, so labels_meta keeps
    // its bytes with the knob off.
    if (p.extra_smear_tran_collection > 0.0) cfg["extra_smear_tran_collection"] = p.extra_smear_tran_collection;
    if (p.extra_smear_tran_induction > 0.0) cfg["extra_smear_tran_induction"] = p.extra_smear_tran_induction;
    // Record which depos fed each blob (per-point refinement in
    // TruthPointLabeler).  Bookkeeping only; off by default.
    cfg["keep_depos"] = p.keep_depos;
    // Fold a negative track key (larg4's "-ancestor" of a dropped EM
    // shower daughter) to the ancestor.  Off = the key stays negative
    // (R1 behaviour); see DepoBlobAssocParams.
    cfg["fold_negative_ids"] = p.fold_negative_ids;
    // Zero-extent depos (backward-drifted from between the response
    // plane and the wires) as points; see DepoBlobAssocParams.  Written
    // only when on: TruthPointLabeler records this configuration in its
    // labels_meta, which keeps its bytes with the knob off.
    if (p.point_zero_width) cfg["point_zero_width"] = true;
    // Truth file (truth/python/dump_truth_table.py) whose sed_ids map a
    // Stage A index-mode depo id to the G4 track id (and sed_energy to
    // the deposit energy).  "" = the depo ids are track ids already.
    cfg["truth_table"] = "";
    return cfg;
}

void Truth::assoc_configure(DepoBlobAssocParams& p, const Configuration& cfg)
{
    p.time_offset = get(cfg, "time_offset", p.time_offset);
    p.nsigma = get(cfg, "nsigma", p.nsigma);
    p.speed = get(cfg, "speed", p.speed);
    p.pindex = get(cfg, "pindex", p.pindex);
    p.id_from_prior = get(cfg, "id_from_prior", p.id_from_prior);
    p.extra_smear_long = get(cfg, "extra_smear_long", p.extra_smear_long);
    p.extra_smear_tran = get(cfg, "extra_smear_tran", p.extra_smear_tran);
    p.extra_smear_tran_collection = get(cfg, "extra_smear_tran_collection", p.extra_smear_tran_collection);
    p.extra_smear_tran_induction = get(cfg, "extra_smear_tran_induction", p.extra_smear_tran_induction);
    if ((p.extra_smear_tran_collection > 0.0 || p.extra_smear_tran_induction > 0.0) && p.extra_smear_tran > 0.0) {
        THROW(ValueError() << errmsg{"extra_smear_tran and the per-plane transverse widths are exclusive"});
    }
    p.keep_depos = get(cfg, "keep_depos", p.keep_depos);
    p.fold_negative_ids = get(cfg, "fold_negative_ids", p.fold_negative_ids);
    p.point_zero_width = get(cfg, "point_zero_width", p.point_zero_width);
    const std::string tfile = get<std::string>(cfg, "truth_table", "");
    if (!tfile.empty()) {
        auto tt = TruthTable::load(tfile);
        if (!tt.sed_ids.empty()) {
            auto map = std::make_shared<std::vector<int>>();
            map->reserve(tt.sed_ids.size());
            for (const auto& [tid, orig] : tt.sed_ids) map->push_back(tid);
            p.id_map = map;
        }
        if (!tt.sed_energy.empty()) {
            p.energy_map = std::make_shared<std::vector<float>>(tt.sed_energy);
        }
    }
    if (p.speed <= 0.0) {
        THROW(ValueError() << errmsg{"no speed provided"});
    }
}

Truth::BlobDepoLabel::BlobDepoLabel()
    : Aux::Logger("BlobDepoLabel", "truth")
{
}
Truth::BlobDepoLabel::~BlobDepoLabel()
{
}

Configuration Truth::BlobDepoLabel::default_configuration() const
{
    return assoc_default_configuration(m_params);
}

void Truth::BlobDepoLabel::configure(const Configuration& cfg)
{
    assoc_configure(m_params, cfg);
}

bool Truth::BlobDepoLabel::operator()(const input_tuple_type& intup,
                                      output_pointer& out)
{
    out = nullptr;
    auto icluster = std::get<0>(intup);
    auto ideposet = std::get<1>(intup);
    if (!icluster or !ideposet) {
        if (icluster) {
            log->warn("join node got non-EOS ICluster {} on EOS", icluster->ident());
        }
        if (ideposet) {
            log->warn("join node got non-EOS IDepoSet {} on EOS", ideposet->ident());
        }
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    auto res = depo_blob_assoc(icluster, *ideposet->depos(), m_params, log);

    size_t nfilled = 0, nmixed = 0;
    for (const auto& [vd, e] : res.blobs) {
        if (e.by_track.empty()) continue;
        ++nfilled;
        if (e.by_track.size() > 1) ++nmixed;
    }
    log->debug("call={} cluster {} deposet {}: nblobs={} ndepos={} filled={} mixed={}",
               m_count, icluster->ident(), ideposet->ident(),
               res.nblobs, res.ndepos, nfilled, nmixed);

    out = filled_cluster(icluster, res);
    ++m_count;
    return true;
}
