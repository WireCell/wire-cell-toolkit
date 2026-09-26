#include "WireCellFlash/SBNDOpFlashFinder.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>

WIRECELL_FACTORY(SBNDOpFlashFinder, WireCell::Flash::SBNDOpFlashFinder,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using Params = Flash::SBNDOpFlashFinder::Params;

namespace {
    struct Hit {
        int channel;
        double time, pe;
    };

    struct Cand {
        std::vector<int> hits;  // indices into the hit list
        int group{-1};          // claimed group it came from; split parts share it
        double time{0}, time_width{0}, first{0}, last{0}, total_pe{0}, y{0}, z{0}, y_width{0}, z_width{0};
        std::vector<double> pes;
    };

    // PE per OpDet, PE-weighted time, half time span, first/last hit, y/z centroid and width.
    void build(Cand& f, const std::vector<Hit>& hits, int nchan, const std::vector<double>& oy,
               const std::vector<double>& oz)
    {
        f.pes.assign(nchan, 0.0);
        f.first = std::numeric_limits<double>::max();
        f.last = -f.first;
        double tsum = 0, pe = 0, sy = 0, sy2 = 0, sz = 0, sz2 = 0;
        for (int k : f.hits) {
            const auto& h = hits[k];
            f.pes[h.channel] += h.pe;
            tsum += h.pe * h.time;
            pe += h.pe;
            f.first = std::min(f.first, h.time);
            f.last = std::max(f.last, h.time);
            sy += h.pe * oy[h.channel];
            sy2 += h.pe * oy[h.channel] * oy[h.channel];
            sz += h.pe * oz[h.channel];
            sz2 += h.pe * oz[h.channel] * oz[h.channel];
        }
        f.total_pe = pe;
        f.time = tsum / pe;
        f.time_width = 0.5 * (f.last - f.first);
        f.y = sy / pe;
        f.z = sz / pe;
        auto width = [&](double s, double s2) {
            const double d = s2 * pe - s * s;
            return d < 0 ? 0.0 : std::sqrt(d) / pe;
        };
        f.y_width = width(sy, sy2);
        f.z_width = width(sz, sz2);
    }

    // Steps 1-2: two offset sets of bins, then claiming biggest candidate first (the larana
    // OpFlashAlg accumulator + AssignHitsToFlash, same order and tie-breaking as OpFlashFinder).
    std::vector<std::vector<int>> accumulate_and_claim(const std::vector<Hit>& hits, const Params& p)
    {
        std::vector<std::vector<int>> groups;
        if (hits.empty()) return groups;
        double tmin = std::numeric_limits<double>::max();
        for (const auto& h : hits) tmin = std::min(tmin, h.time);
        std::vector<double> sum[2];
        std::vector<std::vector<int>> members[2];
        std::vector<int> candidates[2];
        for (size_t k = 0; k < hits.size(); ++k) {
            for (int s = 0; s < 2; ++s) {
                const size_t b = size_t((hits[k].time - tmin + s * p.bin_width / 2.0) / p.bin_width);
                if (b >= sum[s].size()) {
                    sum[s].resize(b + 1, 0.0);
                    members[s].resize(b + 1);
                }
                members[s][b].push_back(k);
                sum[s][b] += hits[k].pe;
                if (sum[s][b] >= p.flash_threshold && sum[s][b] - hits[k].pe < p.flash_threshold)
                    candidates[s].push_back(b);
            }
        }
        std::map<double, std::map<int, std::vector<int>>, std::greater<double>> by_size;
        for (int s = 0; s < 2; ++s)
            for (int b : candidates[s]) by_size[sum[s][b]][s].push_back(b);
        std::vector<char> claimed(hits.size(), 0);
        for (const auto& [pe_unused, per_set] : by_size) {
            for (const auto& [s, bins] : per_set) {
                for (int b : bins) {
                    std::vector<int> mine;
                    double pe = 0;
                    for (int k : members[s][b]) {
                        if (claimed[k]) continue;
                        mine.push_back(k);
                        pe += hits[k].pe;
                    }
                    if (pe < p.flash_threshold) continue;
                    for (int k : mine) claimed[k] = 1;
                    groups.push_back(std::move(mine));
                }
            }
        }
        return groups;
    }

    // Step 3 helper: first place where the hits hold a second prompt burst of light.  Hits in
    // split_bin_ns bins; a spike bin holds >= max(split_spike_pe, split_spike_frac * PE); the
    // first spike is the onset.  A later spike >= split_min_gap_ns after it starts a second
    // flash when its split_burst_ns burst (without its brightest OpDet if split_drop_brightest:
    // a single-PMT artefact is not a flash) is >= max(split_min_pe, split_min_ratio * the onset
    // burst), the split_dip_ns before it hold <= split_dip_frac of that burst, and both parts
    // keep >= flash_threshold.  Returns the cut time (hits at or after it go to the later part).
    bool find_split(const std::vector<int>& fh, const std::vector<Hit>& hits, const Params& p, double& cut)
    {
        if (fh.empty()) return false;
        double t_lo = std::numeric_limits<double>::max(), t_hi = -t_lo, total = 0;
        for (int k : fh) {
            t_lo = std::min(t_lo, hits[k].time);
            t_hi = std::max(t_hi, hits[k].time);
            total += hits[k].pe;
        }
        const double bin = p.split_bin_ns;
        const long nb = long((t_hi - t_lo) / bin) + 1;
        std::vector<double> pe(nb, 0.0);
        for (int k : fh) pe[std::min(nb - 1, long((hits[k].time - t_lo) / bin))] += hits[k].pe;
        auto sum = [&](long a, long b) {
            double s = 0;
            for (long i = std::max(0L, a); i < std::min(nb, b); ++i) s += pe[i];
            return s;
        };
        auto burst = [&](long a, long b) {
            const double s = sum(a, b);
            if (!p.split_drop_brightest) return s;
            const double lo = t_lo + a * bin, hi = t_lo + b * bin;
            std::map<int, double> per;
            for (int k : fh)
                if (hits[k].time >= lo && hits[k].time < hi) per[hits[k].channel] += hits[k].pe;
            double top = 0;
            for (const auto& [ch, v] : per) top = std::max(top, v);
            return s - top;
        };
        const long nburst = std::max(1L, std::lround(p.split_burst_ns / bin));
        const long ngap = std::lround(std::ceil(p.split_min_gap_ns / bin));
        const long ndip = std::lround(p.split_dip_ns / bin);
        const double thr = std::max(p.split_spike_pe, p.split_spike_frac * total);
        long i0 = -1;
        for (long i = 0; i < nb; ++i)
            if (pe[i] >= thr) { i0 = i; break; }
        if (i0 < 0) return false;
        const double p0 = burst(i0, i0 + nburst);
        for (long i = i0 + ngap; i < nb; ++i) {
            if (pe[i] < thr) continue;
            const double b = burst(i, i + nburst);
            if (b < std::max(p.split_min_pe, p.split_min_ratio * p0)) continue;
            if (sum(std::max(i - ndip, i0 + nburst), i) > p.split_dip_frac * b) continue;
            const double early = sum(0, i);
            if (early < p.flash_threshold || total - early < p.flash_threshold) continue;
            cut = t_lo + i * bin;
            return true;
        }
        return false;
    }

    // Join helper: does f contain a sharp burst of its own -- a split_burst_ns window (without
    // its brightest OpDet if split_drop_brightest) holding >= split_min_pe, with the split_dip_ns
    // before it (hits of the earlier flash and f) holding <= split_dip_frac of it?  A tail piece
    // has none; a real new flash does.  Fixed thresholds, not scaled by the flashes' total PE.
    bool has_burst(const std::vector<int>& before, const std::vector<int>& fh, const std::vector<Hit>& hits,
                   const Params& p)
    {
        if (fh.empty()) return false;
        double f0 = std::numeric_limits<double>::max(), f1 = -f0;
        for (int k : fh) {
            f0 = std::min(f0, hits[k].time);
            f1 = std::max(f1, hits[k].time);
        }
        const double bin = p.split_bin_ns, t_lo = f0 - p.split_dip_ns;
        const long nb = long((f1 - t_lo) / bin) + 1;
        std::vector<double> pe(nb, 0.0);
        for (const auto* v : {&before, &fh})
            for (int k : *v) {
                const double t = hits[k].time;
                if (t >= t_lo && t <= f1) pe[std::min(nb - 1, long((t - t_lo) / bin))] += hits[k].pe;
            }
        const long nburst = std::max(1L, std::lround(p.split_burst_ns / bin));
        const long ndip = std::lround(p.split_dip_ns / bin);
        auto sum = [&](long a, long b) {
            double s = 0;
            for (long i = std::max(0L, a); i < std::min(nb, b); ++i) s += pe[i];
            return s;
        };
        for (long i = ndip; i < nb; ++i) {  // burst windows starting inside f
            double b = sum(i, i + nburst);
            if (b < p.split_min_pe) continue;
            if (p.split_drop_brightest) {
                const double lo = t_lo + i * bin, hi = t_lo + (i + nburst) * bin;
                std::map<int, double> per;
                for (int k : fh)
                    if (hits[k].time >= lo && hits[k].time < hi) per[hits[k].channel] += hits[k].pe;
                double top = 0;
                for (const auto& [ch, v] : per) top = std::max(top, v);
                b -= top;
                if (b < p.split_min_pe) continue;
            }
            if (sum(i - ndip, i) <= p.split_dip_frac * b) return true;
        }
        return false;
    }

    std::set<int> lit(const Cand& f, double fired_pe)
    {
        std::set<int> s;
        for (size_t od = 0; od < f.pes.size(); ++od)
            if (f.pes[od] >= fired_pe) s.insert(od);
        return s;
    }

    // SBND flash time (SimpleFlashAlgo candidate bin + FlashT0SelectedChannels).
    double prompt_time(const std::vector<int>& fh, const std::vector<Hit>& hits, const Params& p)
    {
        double t_lo = std::numeric_limits<double>::max();
        for (int k : fh) t_lo = std::min(t_lo, hits[k].time);
        std::map<long, double> bins;
        for (int k : fh) bins[long((hits[k].time - t_lo) / p.prompt_bin_ns)] += hits[k].pe;
        long best = 0;
        double best_pe = -1;
        for (const auto& [b, pe] : bins)  // ascending: ties go to the earliest bin
            if (pe > best_pe) { best = b; best_pe = pe; }
        const double tb = t_lo + best * p.prompt_bin_ns;
        std::vector<std::pair<double, double>> sel;  // (pe, time)
        double pe_sum = 0;
        for (int k : fh) {
            const auto& h = hits[k];
            if (h.time < tb + p.prompt_post_ns && h.time > tb - p.prompt_pre_ns && h.pe > p.prompt_min_hit_pe) {
                sel.emplace_back(h.pe, h.time);
                pe_sum += h.pe;
            }
        }
        if (pe_sum <= 0) return tb;
        std::sort(sel.begin(), sel.end(), std::greater<std::pair<double, double>>());
        double tsum = 0, pe_count = 0;
        int n = 0;
        for (const auto& [pe, t] : sel) {
            pe_count += pe;
            tsum += t;
            ++n;
            if (pe_count / pe_sum > p.prompt_pe_fraction) break;
        }
        return tsum / n;
    }
}

Flash::SBNDOpFlashFinder::SBNDOpFlashFinder()
    : Aux::Logger("SBNDOpFlashFinder", "flash")
{
}

Flash::SBNDOpFlashFinder::~SBNDOpFlashFinder() {}

WireCell::Configuration Flash::SBNDOpFlashFinder::default_configuration() const
{
    const auto& p = m_par;
    Configuration cfg;
    cfg["nchan"] = m_nchan;
    cfg["geom_file"] = m_geom_file;
    cfg["bin_width"] = p.bin_width;
    cfg["flash_threshold"] = p.flash_threshold;
    cfg["pulse_split"] = p.pulse_split;
    cfg["split_bin_ns"] = p.split_bin_ns;
    cfg["split_burst_ns"] = p.split_burst_ns;
    cfg["split_min_gap_us"] = p.split_min_gap_ns / 1000.0;
    cfg["split_spike_pe"] = p.split_spike_pe;
    cfg["split_spike_frac"] = p.split_spike_frac;
    cfg["split_min_ratio"] = p.split_min_ratio;
    cfg["split_min_pe"] = p.split_min_pe;
    cfg["split_dip_ns"] = p.split_dip_ns;
    cfg["split_dip_frac"] = p.split_dip_frac;
    cfg["split_drop_brightest"] = p.split_drop_brightest;
    cfg["join"] = p.join;
    cfg["join_max_gap_us"] = p.join_max_gap_ns / 1000.0;
    cfg["join_pe_ratio"] = p.join_pe_ratio;
    cfg["join_fired_pe"] = p.join_fired_pe;
    cfg["join_lit_frac"] = p.join_lit_frac;
    cfg["join_glow"] = p.join_glow;
    cfg["remove_late_light"] = p.remove_late_light;
    cfg["late_tau_us"] = p.late_tau_ns / 1000.0;
    cfg["late_nsigma"] = p.late_nsigma;
    cfg["min_fired_pds"] = p.min_fired_pds;
    cfg["min_fired_pe"] = p.min_fired_pe;
    cfg["min_total_pe"] = p.min_total_pe;
    cfg["prompt_bin_ns"] = p.prompt_bin_ns;
    cfg["prompt_pre_ns"] = p.prompt_pre_ns;
    cfg["prompt_post_ns"] = p.prompt_post_ns;
    cfg["prompt_min_hit_pe"] = p.prompt_min_hit_pe;
    cfg["prompt_pe_fraction"] = p.prompt_pe_fraction;
    cfg["offset_us"] = m_offset_us;
    cfg["metadata_extra"] = m_metadata_extra;
    return cfg;
}

void Flash::SBNDOpFlashFinder::configure(const WireCell::Configuration& cfg)
{
    auto& p = m_par;
    m_nchan = get(cfg, "nchan", m_nchan);
    m_geom_file = get(cfg, "geom_file", m_geom_file);
    p.bin_width = get(cfg, "bin_width", p.bin_width);
    p.flash_threshold = get(cfg, "flash_threshold", p.flash_threshold);
    p.pulse_split = get(cfg, "pulse_split", p.pulse_split);
    p.split_bin_ns = get(cfg, "split_bin_ns", p.split_bin_ns);
    p.split_burst_ns = get(cfg, "split_burst_ns", p.split_burst_ns);
    p.split_min_gap_ns = 1000.0 * get(cfg, "split_min_gap_us", p.split_min_gap_ns / 1000.0);
    p.split_spike_pe = get(cfg, "split_spike_pe", p.split_spike_pe);
    p.split_spike_frac = get(cfg, "split_spike_frac", p.split_spike_frac);
    p.split_min_ratio = get(cfg, "split_min_ratio", p.split_min_ratio);
    p.split_min_pe = get(cfg, "split_min_pe", p.split_min_pe);
    p.split_dip_ns = get(cfg, "split_dip_ns", p.split_dip_ns);
    p.split_dip_frac = get(cfg, "split_dip_frac", p.split_dip_frac);
    p.split_drop_brightest = get(cfg, "split_drop_brightest", p.split_drop_brightest);
    p.join = get(cfg, "join", p.join);
    p.join_max_gap_ns = 1000.0 * get(cfg, "join_max_gap_us", p.join_max_gap_ns / 1000.0);
    p.join_pe_ratio = get(cfg, "join_pe_ratio", p.join_pe_ratio);
    p.join_fired_pe = get(cfg, "join_fired_pe", p.join_fired_pe);
    p.join_lit_frac = get(cfg, "join_lit_frac", p.join_lit_frac);
    p.join_glow = get(cfg, "join_glow", p.join_glow);
    if (p.join_glow != "any" && p.join_glow != "only" && p.join_glow != "never") {
        raise<ValueError>("SBNDOpFlashFinder: join_glow must be any, only or never");
    }
    p.remove_late_light = get(cfg, "remove_late_light", p.remove_late_light);
    p.late_tau_ns = 1000.0 * get(cfg, "late_tau_us", p.late_tau_ns / 1000.0);
    p.late_nsigma = get(cfg, "late_nsigma", p.late_nsigma);
    p.min_fired_pds = get(cfg, "min_fired_pds", p.min_fired_pds);
    p.min_fired_pe = get(cfg, "min_fired_pe", p.min_fired_pe);
    p.min_total_pe = get(cfg, "min_total_pe", p.min_total_pe);
    p.prompt_bin_ns = get(cfg, "prompt_bin_ns", p.prompt_bin_ns);
    p.prompt_pre_ns = get(cfg, "prompt_pre_ns", p.prompt_pre_ns);
    p.prompt_post_ns = get(cfg, "prompt_post_ns", p.prompt_post_ns);
    p.prompt_min_hit_pe = get(cfg, "prompt_min_hit_pe", p.prompt_min_hit_pe);
    p.prompt_pe_fraction = get(cfg, "prompt_pe_fraction", p.prompt_pe_fraction);
    m_offset_us = get(cfg, "offset_us", m_offset_us);
    if (cfg.isMember("metadata_extra")) m_metadata_extra = cfg["metadata_extra"];
    if (p.bin_width <= 0 || p.split_bin_ns <= 0 || p.prompt_bin_ns <= 0) {
        raise<ValueError>("SBNDOpFlashFinder: bin widths must be > 0");
    }

    m_opdet_y.assign(m_nchan, 0.0);
    m_opdet_z.assign(m_nchan, 0.0);
    if (m_geom_file.empty()) {
        log->warn("no geom_file: flash y/z centroids will be 0");
        return;
    }
    auto jgeom = Persist::load(m_geom_file);
    for (const auto& jod : jgeom["opdets"]) {
        const int od = jod["opdet"].asInt();
        if (od < 0 || od >= m_nchan) continue;
        m_opdet_y[od] = jod["y"].asDouble();
        m_opdet_z[od] = jod["z"].asDouble();
    }
}

bool Flash::SBNDOpFlashFinder::operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    ++m_count;
    const auto& p = m_par;

    ITensor::pointer hits_ten = nullptr;
    for (const auto& ten : *in->tensors()) {
        if (ten->metadata()["name"].asString() == "ophits") {
            hits_ten = ten;
            break;
        }
    }
    if (!hits_ten) hits_ten = in->tensors()->at(0);
    const size_t nrow = hits_ten->shape()[0];
    const size_t ncol = hits_ten->shape()[1];
    if (ncol < 9) raise<ValueError>("SBNDOpFlashFinder: ophits tensor has %d < 9 columns", (int) ncol);
    const double* H = (const double*) hits_ten->data();

    // Hits on channels outside [0, nchan) are ignored (kept in the output with no flash).
    std::vector<Hit> hits;
    std::vector<size_t> row_of;  // hit index -> input row
    for (size_t r = 0; r < nrow; ++r) {
        const double* row = H + r * ncol;
        const int ch = int(row[0]);
        if (ch < 0 || ch >= m_nchan) continue;
        hits.push_back(Hit{ch, row[1], row[5]});
        row_of.push_back(r);
    }

    auto make = [&](std::vector<int> fh, int group) {
        Cand f;
        f.hits = std::move(fh);
        f.group = group;
        build(f, hits, m_nchan, m_opdet_y, m_opdet_z);
        return f;
    };
    auto by_first = [](const Cand& a, const Cand& b) { return a.first < b.first; };

    // 1-2. accumulate + claim
    auto groups = accumulate_and_claim(hits, p);
    std::vector<Cand> flashes;
    for (size_t g = 0; g < groups.size(); ++g) flashes.push_back(make(std::move(groups[g]), g));

    // 3. pulse split (parts keep the group id of the flash they came from)
    if (p.pulse_split) {
        std::vector<Cand> parts;
        for (auto& f : flashes) {
            std::vector<int> rest = f.hits;
            double cut = 0;
            bool split = false;
            while (find_split(rest, hits, p, cut)) {
                std::vector<int> early, late;
                for (int k : rest) (hits[k].time < cut ? early : late).push_back(k);
                parts.push_back(make(std::move(early), f.group));
                rest = std::move(late);
                split = true;
            }
            if (split) parts.push_back(make(std::move(rest), f.group));
            else parts.push_back(std::move(f));
        }
        flashes = std::move(parts);
    }

    // 4. join each flash to the nearest earlier one (see header)
    std::sort(flashes.begin(), flashes.end(), by_first);
    if (p.join && flashes.size() > 1) {
        std::vector<Cand> kept;
        for (auto& f : flashes) {
            bool joined = false;
            if (!kept.empty()) {
                Cand& i = kept.back();
                // most of f's light on OpDets the earlier flash lights (not all: the two flashes'
                // tails overlap, so a piece carries a few hits of other flashes)
                const auto lit_i = lit(i, p.join_fired_pe);
                double on_lit = 0;
                for (int od : lit_i) on_lit += f.pes[od];
                const bool subset = f.total_pe > 0 && on_lit >= p.join_lit_frac * f.total_pe;
                const bool pe_ok = p.join_pe_ratio <= 0 || f.total_pe <= p.join_pe_ratio * i.total_pe;
                // the earlier flash's slow tail must explain f: expected = i's PE times the share
                // of an exponential (late_tau) starting at i's onset that falls in f's time span,
                // normalised to the span i covers.  Lenient on purpose (i's fast light counted as
                // tail); a piece with significantly more light is a separate flash.
                // join_glow: "any" = no glow test; "only" = join only pieces the earlier flash's
                // glow explains; "never" = join only pieces with more light than that glow
                // (glow pieces are then left to step 5)
                bool tail_ok = true;
                if (p.join_glow != "any") {
                    const double t0 = prompt_time(i.hits, hits, p);
                    const double tau = p.late_tau_ns;
                    const double norm = 1.0 - std::exp(-(i.last - t0) / tau);
                    const double hyp = norm > 0 ? i.total_pe *
                        (std::exp(-(f.first - t0) / tau) - std::exp(-(f.last - t0) / tau)) / norm : 0.0;
                    const bool glow = hyp > 0 && (f.total_pe - hyp) / std::sqrt(hyp) < p.late_nsigma;
                    tail_ok = (p.join_glow == "only") ? glow : !glow;
                }
                if (i.group != f.group && f.first - i.last <= p.join_max_gap_ns && pe_ok && subset &&
                    tail_ok) {
                    // a real new flash has a sharp burst of its own; a tail piece does not
                    const bool new_flash = has_burst(i.hits, f.hits, hits, p);
                    if (!new_flash) {
                        std::vector<int> both = i.hits;
                        both.insert(both.end(), f.hits.begin(), f.hits.end());
                        i = make(std::move(both), i.group);
                        joined = true;
                    }
                }
            }
            if (!joined) kept.push_back(std::move(f));
        }
        flashes = std::move(kept);
    }

    // 5. remove late light (larana RemoveLateLight; PE-weighted times at this stage)
    if (p.remove_late_light && flashes.size() > 1) {
        std::sort(flashes.begin(), flashes.end(),
                  [](const Cand& a, const Cand& b) { return a.time < b.time; });
        std::vector<char> remove(flashes.size(), 0);
        for (size_t i = 0; i < flashes.size(); ++i) {
            for (size_t j = i + 1; j < flashes.size(); ++j) {
                if (remove[j]) continue;
                const auto& fi = flashes[i];
                const auto& fj = flashes[j];
                if (fi.time > fj.time || fi.time_width <= 0) continue;
                // never between two parts of one split (the split found a new burst there), and a
                // flash cannot be the tail of a dimmer one
                if (fi.group == fj.group || fj.total_pe >= fi.total_pe) continue;
                const double hyp = fi.total_pe * fj.time_width / fi.time_width *
                                   std::exp(-(fj.time - fi.time) / p.late_tau_ns);
                if (hyp <= 0) continue;
                if ((fj.total_pe - hyp) / std::sqrt(hyp) < p.late_nsigma) remove[j] = 1;
            }
        }
        std::vector<Cand> kept;
        for (size_t i = 0; i < flashes.size(); ++i)
            if (!remove[i]) kept.push_back(std::move(flashes[i]));
        flashes = std::move(kept);
    }

    // 6. quality cut
    {
        std::vector<Cand> kept;
        for (auto& f : flashes) {
            int npd = 0;
            for (double v : f.pes) npd += (v >= p.min_fired_pe);
            if (npd >= p.min_fired_pds && f.total_pe >= p.min_total_pe) kept.push_back(std::move(f));
        }
        flashes = std::move(kept);
    }

    // 7. SBND prompt time, then time order
    for (auto& f : flashes) f.time = prompt_time(f.hits, hits, p);
    std::stable_sort(flashes.begin(), flashes.end(),
                     [](const Cand& a, const Cand& b) { return a.time < b.time; });

    // Output: the OpFlashFinder tensor-set schema.
    const size_t nflash = flashes.size();
    const size_t mcol = 1 + m_nchan;
    std::vector<double> matrix(nflash * mcol, 0.0), summary(nflash * 8, 0.0), ohits(nrow * 9);
    for (size_t r = 0; r < nrow; ++r) {
        std::copy(H + r * ncol, H + r * ncol + 9, &ohits[r * 9]);
        ohits[r * 9 + 7] = -1;
    }
    for (size_t f = 0; f < nflash; ++f) {
        const auto& fl = flashes[f];
        matrix[f * mcol] = fl.time;
        std::copy(fl.pes.begin(), fl.pes.end(), &matrix[f * mcol + 1]);
        double* s = &summary[f * 8];
        s[0] = f;
        s[1] = fl.total_pe;
        s[2] = fl.y;
        s[3] = fl.z;
        s[4] = fl.y_width;
        s[5] = fl.z_width;
        s[6] = -1;
        s[7] = fl.hits.size();
        for (int k : fl.hits) ohits[row_of[k] * 9 + 7] = f;
    }
    ITensor::vector* tensors = new ITensor::vector;
    auto add = [&](const char* name, size_t nr, size_t nc, const std::vector<double>& v) {
        Configuration md;
        md["name"] = name;
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{nr, nc}, v.data(), md));
    };
    add("opflash", nflash, mcol, matrix);
    add("flash_summary", nflash, 8, summary);
    add("ophits", nrow, 9, ohits);

    Configuration md = in->metadata();
    md["producer"] = "wct-flash-sbnd";
    md["nchan"] = m_nchan;
    md["offset_us"] = m_offset_us;
    if (!m_metadata_extra.isNull())
        for (const auto& key : m_metadata_extra.getMemberNames()) md[key] = m_metadata_extra[key];
    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), md, ITensor::shared_vector(tensors));
    log->debug("set {}: {} flashes from {} hits", in->ident(), nflash, hits.size());
    return true;
}
