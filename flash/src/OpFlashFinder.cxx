#include "WireCellFlash/OpFlashFinder.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>

WIRECELL_FACTORY(OpFlashFinder, WireCell::Flash::OpFlashFinder,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;

namespace {
    // One ophit row (see OpHitFinder output schema).
    struct Hit {
        int channel;
        double peak_time, width, area, amplitude, pe, start_time, fast_to_total;
    };

    // Port of larana opdet::RunFlashFinder and helpers, on Hit rows.
    // Times in WCT ns.

    void fill_accumulator(unsigned index, unsigned hit_index, double pe, double flash_threshold,
                          std::vector<double>& binned, std::vector<std::vector<int>>& contributors,
                          std::vector<int>& flashes_in_accumulator)
    {
        contributors.at(index).push_back(hit_index);
        binned.at(index) += pe;
        if (binned.at(index) >= flash_threshold && (binned.at(index) - pe) < flash_threshold) {
            flashes_in_accumulator.push_back(index);
        }
    }

    void assign_hits_to_flash(const std::vector<int>& flashes1, const std::vector<int>& flashes2,
                              const std::vector<double>& binned1, const std::vector<double>& binned2,
                              const std::vector<std::vector<int>>& contributors1,
                              const std::vector<std::vector<int>>& contributors2,
                              const std::vector<Hit>& hits,
                              std::vector<std::vector<int>>& hits_per_flash,
                              double flash_threshold)
    {
        // FlashesBySize[PE][accumulator] = bins, walked largest first.
        std::map<double, std::map<int, std::vector<int>>, std::greater<double>> flashes_by_size;
        for (int bin : flashes1) flashes_by_size[binned1.at(bin)][1].push_back(bin);
        for (int bin : flashes2) flashes_by_size[binned2.at(bin)][2].push_back(bin);

        std::vector<int> claimed(hits.size(), -1);
        for (const auto& it_flash : flashes_by_size) {
            for (const auto& it_acc : it_flash.second) {
                const auto& contributors = (it_acc.first == 1) ? contributors1 : contributors2;
                for (int bin : it_acc.second) {
                    std::vector<int> hits_this_flash;
                    for (int hit_index : contributors.at(bin)) {
                        if (claimed.at(hit_index) == -1) hits_this_flash.push_back(hit_index);
                    }
                    double pe = 0;
                    for (int hit_index : hits_this_flash) pe += hits.at(hit_index).pe;
                    if (pe < flash_threshold) continue;
                    hits_per_flash.push_back(hits_this_flash);
                    for (int hit_index : hits_this_flash) {
                        claimed.at(hit_index) = hits_per_flash.size() - 1;
                    }
                }
            }
        }
    }

    void refine_hits_in_flash(const std::vector<int>& hits_this_flash, const std::vector<Hit>& hits,
                              std::vector<std::vector<int>>& refined, double width_tolerance,
                              double flash_threshold)
    {
        std::map<double, std::vector<int>, std::greater<double>> hits_by_size;
        for (int hit_index : hits_this_flash) hits_by_size[hits.at(hit_index).pe].push_back(hit_index);

        std::vector<bool> used(hits.size(), false);
        while (true) {
            std::vector<int> flash_hits;
            double pe = 0, tmax = 0, tmin = 0;
            // Seed with the biggest unused hit.
            for (const auto& it : hits_by_size) {
                for (int hit_index : it.second) {
                    if (used.at(hit_index)) continue;
                    const auto& h = hits.at(hit_index);
                    pe = h.pe;
                    tmax = h.peak_time + 0.5 * h.width;
                    tmin = h.peak_time - 0.5 * h.width;
                    flash_hits.push_back(hit_index);
                    used.at(hit_index) = true;
                    goto seeded;
                }
            }
          seeded:
            if (flash_hits.empty()) return;

            // Collect hits within the width tolerance until stable.
            size_t nlast = 0;
            while (nlast < flash_hits.size()) {
                nlast = flash_hits.size();
                for (const auto& it : hits_by_size) {
                    for (int hit_index : it.second) {
                        if (used.at(hit_index)) continue;
                        const auto& h = hits.at(hit_index);
                        const double hw = 0.5 * h.width;
                        const double ft = 0.5 * (tmax + tmin);
                        const double fw = 0.5 * (tmax - tmin);
                        if (std::abs(h.peak_time - ft) > width_tolerance * (hw + fw)) continue;
                        flash_hits.push_back(hit_index);
                        tmax = std::max(tmax, h.peak_time + hw);
                        tmin = std::min(tmin, h.peak_time - hw);
                        pe += h.pe;
                        used.at(hit_index) = true;
                    }
                }
            }

            if (pe >= flash_threshold) {
                refined.push_back(flash_hits);
            }
            else if (flash_hits.size() > 1) {
                // Release all but the seed for possible reuse.
                for (size_t i = 1; i < flash_hits.size(); ++i) used.at(flash_hits[i]) = false;
            }
        }
    }

    struct FlashSummary {
        double time, time_width, total_pe, y, y_width, z, z_width;
        std::vector<double> pes;
    };

    FlashSummary construct_flash(const std::vector<int>& flash_hits, const std::vector<Hit>& hits,
                                 int nchan, const std::vector<double>& opdet_y,
                                 const std::vector<double>& opdet_z)
    {
        FlashSummary fs{};
        fs.pes.assign(nchan, 0.0);
        double max_time = -std::numeric_limits<double>::max();
        double min_time = std::numeric_limits<double>::max();
        double sumy = 0, sumy2 = 0, sumz = 0, sumz2 = 0;
        for (int hit_index : flash_hits) {
            const auto& h = hits.at(hit_index);
            max_time = std::max(max_time, h.peak_time);
            min_time = std::min(min_time, h.peak_time);
            fs.time += h.pe * h.peak_time;
            fs.total_pe += h.pe;
            if (h.channel >= 0 && h.channel < nchan) fs.pes[h.channel] += h.pe;
            const double y = opdet_y.at(h.channel);
            const double z = opdet_z.at(h.channel);
            sumy += h.pe * y;
            sumy2 += h.pe * y * y;
            sumz += h.pe * z;
            sumz2 += h.pe * z * z;
        }
        fs.time /= fs.total_pe;
        fs.time_width = 0.5 * (max_time - min_time);
        fs.y = sumy / fs.total_pe;
        fs.z = sumz / fs.total_pe;
        auto width = [&](double sum, double sum2) {
            const double d = sum2 * fs.total_pe - sum * sum;
            return (d < 0) ? 0.0 : std::sqrt(d) / fs.total_pe;
        };
        fs.y_width = width(sumy, sumy2);
        fs.z_width = width(sumz, sumz2);
        return fs;
    }

    void remove_late_light(std::vector<FlashSummary>& flashes,
                           std::vector<std::vector<int>>& hits_per_flash)
    {
        // Sort flashes (and their hit lists) by time, then drop any
        // flash consistent with late light of an earlier one.
        std::vector<int> order(flashes.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return flashes[a].time < flashes[b].time; });
        std::vector<FlashSummary> sorted_flashes;
        std::vector<std::vector<int>> sorted_hits;
        for (int i : order) {
            sorted_flashes.push_back(std::move(flashes[i]));
            sorted_hits.push_back(std::move(hits_per_flash[i]));
        }
        flashes = std::move(sorted_flashes);
        hits_per_flash = std::move(sorted_hits);

        std::vector<bool> remove(flashes.size(), false);
        const double argon_tau = 1600.0;  // ns ("1.6" in the us-based original)
        for (size_t i = 0; i < flashes.size(); ++i) {
            for (size_t j = i + 1; j < flashes.size(); ++j) {
                if (remove[j]) continue;
                const auto& fi = flashes[i];
                const auto& fj = flashes[j];
                if (fi.time > fj.time) continue;  // likelihood 1e6, keep
                const double hyp_pe = fi.total_pe * fj.time_width / fi.time_width *
                                      std::exp(-(fj.time - fi.time) / argon_tau);
                const double nsigma = (fj.total_pe - hyp_pe) / std::sqrt(hyp_pe);
                if (nsigma < 3.0) remove[j] = true;
            }
        }
        for (int i = remove.size() - 1; i >= 0; --i) {
            if (remove[i]) {
                flashes.erase(flashes.begin() + i);
                hits_per_flash.erase(hits_per_flash.begin() + i);
            }
        }
    }

    // Knobs + precomputed per-OpDet grid for the optional refinement pass.
    struct RefineParams {
        bool   enabled   = false;
        double window_ns = 8000.0;   // = m_refine_window_us * 1000 (times are ns)
        double pe_ratio  = 0.5;      // later_pe <= pe_ratio * earlier_pe
        int    max_fired = 2;        // later flash fired-PD count cap
        double fired_pe  = 0.5;      // pes[od] >= fired_pe counts as "fired"
        bool   subset_merge = false; // bypass max_fired when lit_j subset of lit_i
        // [nchan] grid coords within each cathode side (-1 if unmapped).
        const std::vector<int>* row  = nullptr;
        const std::vector<int>* col  = nullptr;
        const std::vector<int>* side = nullptr;
    };

    // Optional "flash refinement": merge a later, dim, few-PD flash into an
    // earlier one whose lit OpDets are physically adjacent.  Walks flashes in
    // time order and, for each earlier flash i, absorbs any later flash j that
    // is (1) within the time window, (2) dim relative to i, (3) lit on only a
    // few OpDets (or, with subset_merge, lit on only OpDets that i already
    // lights -- a tail of a bright extended parent), and (4) every lit OpDet of
    // j is the same as or an 8-neighbour
    // (Chebyshev<=1, same side) of a lit OpDet of i.  Each merge recomputes i
    // via construct_flash, so it cascades: a third flash is tested against the
    // already-grown i.  flashes and refined are kept parallel (same discipline
    // as remove_late_light).  No-op when !rp.enabled => bit-identical default.
    void refine_flashes(std::vector<FlashSummary>& flashes,
                        std::vector<std::vector<int>>& refined,
                        const std::vector<Hit>& hits, int nchan,
                        const std::vector<double>& opdet_y,
                        const std::vector<double>& opdet_z,
                        const RefineParams& rp)
    {
        if (!rp.enabled || flashes.size() < 2) return;

        // Sort flashes (and their hit lists) by time -- the cascade and the
        // window break both need time order, and this must hold even when
        // remove_late_light was disabled (then flashes arrive by size).
        {
            std::vector<int> order(flashes.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&](int a, int b) { return flashes[a].time < flashes[b].time; });
            std::vector<FlashSummary> sf;
            std::vector<std::vector<int>> sr;
            for (int k : order) {
                sf.push_back(std::move(flashes[k]));
                sr.push_back(std::move(refined[k]));
            }
            flashes = std::move(sf);
            refined = std::move(sr);
        }

        const auto& row = *rp.row;
        const auto& col = *rp.col;
        const auto& side = *rp.side;

        auto fired_pds = [&](const FlashSummary& f) {
            std::vector<int> v;
            for (int od = 0; od < nchan; ++od)
                if (f.pes[od] >= rp.fired_pe) v.push_back(od);
            return v;
        };
        auto adjacent = [&](int oda, int odb) {  // 8-neighbour, same side
            return side[oda] >= 0 && side[oda] == side[odb] &&
                   std::max(std::abs(row[oda] - row[odb]),
                            std::abs(col[oda] - col[odb])) <= 1;
        };

        for (size_t i = 0; i < flashes.size(); ++i) {
            std::vector<int> fired_i = fired_pds(flashes[i]);
            size_t j = i + 1;
            while (j < flashes.size()) {
                if (flashes[j].time - flashes[i].time > rp.window_ns) break;  // sorted => done

                const auto& fj = flashes[j];
                std::vector<int> fired_j = fired_pds(fj);

                // subset escape: j lights only PDs that i already lights, so it
                // bypasses the few-PD cap (a tail of an extended bright parent).
                bool subset = rp.subset_merge;
                for (int odj : fired_j) {
                    if (!subset) break;
                    if (std::find(fired_i.begin(), fired_i.end(), odj) == fired_i.end())
                        subset = false;
                }
                bool ok = (fj.total_pe <= rp.pe_ratio * flashes[i].total_pe)  // dim
                       && ((int) fired_j.size() >= 1)                         // a real small flash
                       && (((int) fired_j.size() <= rp.max_fired) || subset); // few PDs (or subset)
                if (ok) {
                    for (int odj : fired_j) {  // every lit j adjacent to some lit i
                        bool near = false;
                        for (int odi : fired_i)
                            if (adjacent(odi, odj)) { near = true; break; }
                        if (!near) { ok = false; break; }
                    }
                }

                if (ok) {
                    refined[i].insert(refined[i].end(),
                                      refined[j].begin(), refined[j].end());
                    flashes[i] = construct_flash(refined[i], hits, nchan, opdet_y, opdet_z);
                    fired_i = fired_pds(flashes[i]);  // i grew: refresh
                    flashes.erase(flashes.begin() + j);
                    refined.erase(refined.begin() + j);
                    // do NOT advance j: the next flash slides into index j and
                    // is tested against the grown i (cascade).  i stays fixed
                    // while j advances, so i can absorb a satellite past an
                    // intervening non-merged flash; the spatial+ratio gates keep
                    // that benign.
                }
                else {
                    ++j;
                }
            }
        }
    }

    // Run the full larana flash-finding pipeline over one subset of hits
    // (given as global indices into `hits`): accumulators -> hit claiming
    // -> width refinement -> flash construction -> late-light removal.
    // `out_refined` is returned in the same global hit indices.  Calling
    // this once with all indices reproduces the original all-OpDet result
    // exactly; calling it per cathode side builds optically-independent
    // flashes (see OpFlashFinder::operator()).
    void find_flashes(const std::vector<Hit>& hits, const std::vector<int>& subset,
                      double bin_width, double flash_threshold, double width_tolerance,
                      bool remove_late, int nchan,
                      const std::vector<double>& opdet_y, const std::vector<double>& opdet_z,
                      const RefineParams& rp,
                      std::vector<FlashSummary>& out_flashes,
                      std::vector<std::vector<int>>& out_refined)
    {
        out_flashes.clear();
        out_refined.clear();
        if (subset.empty()) return;

        double min_time = std::numeric_limits<double>::max();
        for (int r : subset) min_time = std::min(min_time, hits[r].peak_time);

        size_t initial = 6400;
        std::vector<double> binned1(initial), binned2(initial);
        std::vector<std::vector<int>> contributors1(initial), contributors2(initial);
        std::vector<int> flashes1, flashes2;
        for (int r : subset) {
            const auto& h = hits[r];
            const unsigned i1 = unsigned((h.peak_time - min_time) / bin_width);
            const unsigned i2 = unsigned((h.peak_time - min_time + bin_width / 2.0) / bin_width);
            if (i2 >= binned1.size()) {
                const size_t n = i2 * 1.2 + 1;
                binned1.resize(n);
                binned2.resize(n);
                contributors1.resize(n);
                contributors2.resize(n);
            }
            fill_accumulator(i1, r, h.pe, flash_threshold, binned1, contributors1, flashes1);
            fill_accumulator(i2, r, h.pe, flash_threshold, binned2, contributors2, flashes2);
        }

        std::vector<std::vector<int>> hits_per_flash;
        assign_hits_to_flash(flashes1, flashes2, binned1, binned2, contributors1, contributors2,
                             hits, hits_per_flash, flash_threshold);

        for (const auto& flash_hits : hits_per_flash) {
            refine_hits_in_flash(flash_hits, hits, out_refined, width_tolerance, flash_threshold);
        }

        for (const auto& flash_hits : out_refined) {
            out_flashes.push_back(construct_flash(flash_hits, hits, nchan, opdet_y, opdet_z));
        }

        if (remove_late) {
            remove_late_light(out_flashes, out_refined);
        }

        // Optional merge of over-split satellite flashes (no-op if disabled).
        refine_flashes(out_flashes, out_refined, hits, nchan, opdet_y, opdet_z, rp);
    }
}

Flash::OpFlashFinder::OpFlashFinder()
    : Aux::Logger("OpFlashFinder", "flash")
{
}

Flash::OpFlashFinder::~OpFlashFinder() {}

WireCell::Configuration Flash::OpFlashFinder::default_configuration() const
{
    Configuration cfg;
    cfg["nchan"] = m_nchan;
    cfg["geom_file"] = m_geom_file;
    cfg["bin_width"] = m_bin_width;
    cfg["flash_threshold"] = m_flash_threshold;
    cfg["width_tolerance"] = m_width_tolerance;
    cfg["remove_late_light"] = m_remove_late_light;
    cfg["group_by_side"] = m_group_by_side;
    cfg["flash_refine"] = m_flash_refine;
    cfg["refine_window_us"] = m_refine_window_us;
    cfg["refine_pe_ratio"] = m_refine_pe_ratio;
    cfg["refine_max_fired"] = m_refine_max_fired;
    cfg["refine_fired_pe"] = m_refine_fired_pe;
    cfg["refine_subset_merge"] = m_refine_subset_merge;
    cfg["offset_us"] = m_offset_us;
    return cfg;
}

void Flash::OpFlashFinder::configure(const WireCell::Configuration& cfg)
{
    m_nchan = get(cfg, "nchan", m_nchan);
    m_geom_file = get(cfg, "geom_file", m_geom_file);
    m_bin_width = get(cfg, "bin_width", m_bin_width);
    m_flash_threshold = get(cfg, "flash_threshold", m_flash_threshold);
    m_width_tolerance = get(cfg, "width_tolerance", m_width_tolerance);
    m_remove_late_light = get(cfg, "remove_late_light", m_remove_late_light);
    m_group_by_side = get(cfg, "group_by_side", m_group_by_side);
    m_flash_refine = get(cfg, "flash_refine", m_flash_refine);
    m_refine_window_us = get(cfg, "refine_window_us", m_refine_window_us);
    m_refine_pe_ratio = get(cfg, "refine_pe_ratio", m_refine_pe_ratio);
    m_refine_max_fired = get(cfg, "refine_max_fired", m_refine_max_fired);
    m_refine_fired_pe = get(cfg, "refine_fired_pe", m_refine_fired_pe);
    m_refine_subset_merge = get(cfg, "refine_subset_merge", m_refine_subset_merge);
    m_offset_us = get(cfg, "offset_us", m_offset_us);

    auto jgeom = Persist::load(m_geom_file);
    m_opdet_x.assign(m_nchan, 0.0);
    m_opdet_y.assign(m_nchan, 0.0);
    m_opdet_z.assign(m_nchan, 0.0);
    for (const auto& jod : jgeom["opdets"]) {
        const int od = jod["opdet"].asInt();
        if (od < 0 or od >= m_nchan) continue;
        m_opdet_x[od] = jod["x"].asDouble();
        m_opdet_y[od] = jod["y"].asDouble();
        m_opdet_z[od] = jod["z"].asDouble();
    }

    // Build the per-OpDet grid (row = y-rank, col = z-rank, side) once, for the
    // 8-neighbour adjacency test in refine_flashes.  Each cathode side is a
    // regular 10(y) x 8(z) grid; rank the distinct y/z per side with a 1 mm
    // tolerance so float noise never splits a rank (real pitch ~600/~500 mm).
    m_opdet_row.assign(m_nchan, -1);
    m_opdet_col.assign(m_nchan, -1);
    m_opdet_side.assign(m_nchan, -1);
    const double eps = 1.0;  // mm
    auto add_level = [&](std::vector<double>& lv, double v) {
        for (double u : lv) if (std::abs(u - v) <= eps) return;
        lv.push_back(v);
    };
    auto rank_of = [&](double v, const std::vector<double>& levels) {
        for (size_t k = 0; k < levels.size(); ++k)
            if (std::abs(v - levels[k]) <= eps) return (int) k;
        return -1;
    };
    for (int s = 0; s < 2; ++s) {  // side 0: x >= 0, side 1: x < 0
        std::vector<int> ods;
        for (int od = 0; od < m_nchan; ++od)
            if ((m_opdet_x[od] >= 0.0) == (s == 0)) ods.push_back(od);
        std::vector<double> ylev, zlev;
        for (int od : ods) { add_level(ylev, m_opdet_y[od]); add_level(zlev, m_opdet_z[od]); }
        std::sort(ylev.begin(), ylev.end());
        std::sort(zlev.begin(), zlev.end());
        for (int od : ods) {
            m_opdet_row[od] = rank_of(m_opdet_y[od], ylev);
            m_opdet_col[od] = rank_of(m_opdet_z[od], zlev);
            m_opdet_side[od] = s;
        }
    }
}

bool Flash::OpFlashFinder::operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    ++m_count;

    // Locate the ophits tensor (by metadata name, else tensor 0).
    ITensor::pointer hits_ten = nullptr;
    for (const auto& ten : *in->tensors()) {
        if (ten->metadata()["name"].asString() == "ophits") {
            hits_ten = ten;
            break;
        }
    }
    if (!hits_ten) hits_ten = in->tensors()->at(0);
    const auto shape = hits_ten->shape();
    const size_t nhit = shape[0];
    const size_t ncol = shape[1];
    const double* H = (const double*) hits_ten->data();

    std::vector<Hit> hits(nhit);
    for (size_t r = 0; r < nhit; ++r) {
        const double* row = H + r * ncol;
        hits[r] = Hit{int(row[0]), row[1], row[2], row[3], row[4], row[5], row[6], row[8]};
    }

    // Build flashes either across all OpDets (larana default) or
    // independently per cathode side (PDHD: the cathode is opaque, so the
    // two drift volumes are optically independent and a flash belongs to
    // exactly one of them).  The all-OpDet path is bit-identical to the
    // original single-pass code.
    RefineParams rp;
    rp.enabled = m_flash_refine;
    rp.window_ns = m_refine_window_us * 1000.0;  // us -> WCT ns (cf. bin_width)
    rp.pe_ratio = m_refine_pe_ratio;
    rp.max_fired = m_refine_max_fired;
    rp.fired_pe = m_refine_fired_pe;
    rp.subset_merge = m_refine_subset_merge;
    rp.row = &m_opdet_row;
    rp.col = &m_opdet_col;
    rp.side = &m_opdet_side;

    std::vector<FlashSummary> flashes;
    std::vector<std::vector<int>> refined;
    if (m_group_by_side) {
        std::vector<int> sideA, sideB;  // A: x >= 0, B: x < 0
        for (size_t r = 0; r < nhit; ++r) {
            const int ch = hits[r].channel;
            const double x = (ch >= 0 && ch < m_nchan) ? m_opdet_x[ch] : 0.0;
            (x >= 0.0 ? sideA : sideB).push_back((int) r);
        }
        for (const std::vector<int>* sub : {&sideA, &sideB}) {
            std::vector<FlashSummary> fs;
            std::vector<std::vector<int>> rf;
            find_flashes(hits, *sub, m_bin_width, m_flash_threshold, m_width_tolerance,
                         m_remove_late_light, m_nchan, m_opdet_y, m_opdet_z, rp, fs, rf);
            for (size_t k = 0; k < fs.size(); ++k) {
                flashes.push_back(std::move(fs[k]));
                refined.push_back(std::move(rf[k]));
            }
        }
        // Keep flash indices monotonic in time across the two sides
        // (no-op when one side is empty, as in the 2024 single-side data).
        std::vector<int> order(flashes.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return flashes[a].time < flashes[b].time; });
        std::vector<FlashSummary> sf;
        std::vector<std::vector<int>> sr;
        for (int i : order) {
            sf.push_back(std::move(flashes[i]));
            sr.push_back(std::move(refined[i]));
        }
        flashes = std::move(sf);
        refined = std::move(sr);
    }
    else {
        std::vector<int> all(nhit);
        std::iota(all.begin(), all.end(), 0);
        find_flashes(hits, all, m_bin_width, m_flash_threshold, m_width_tolerance,
                     m_remove_late_light, m_nchan, m_opdet_y, m_opdet_z, rp, flashes, refined);
    }

    // Emit the opflash tensor-set schema (design.md §3.4).
    const size_t nflash = flashes.size();
    const size_t mcol = 1 + m_nchan;
    std::vector<double> matrix(nflash * mcol, 0.0);
    std::vector<double> summary(nflash * 8, 0.0);
    std::vector<double> ohits(nhit * 9);
    for (size_t r = 0; r < nhit; ++r) {
        std::copy(H + r * ncol, H + r * ncol + 9, &ohits[r * 9]);
    }
    for (size_t f = 0; f < nflash; ++f) {
        const auto& fs = flashes[f];
        matrix[f * mcol] = fs.time;
        for (int c = 0; c < m_nchan; ++c) matrix[f * mcol + 1 + c] = fs.pes[c];
        double* srow = &summary[f * 8];
        srow[0] = f;
        srow[1] = fs.total_pe;
        srow[2] = fs.y;
        srow[3] = fs.z;
        srow[4] = fs.y_width;
        srow[5] = fs.z_width;
        srow[6] = -1;  // absolute DTS time not tracked here
        srow[7] = refined[f].size();
        for (int hit_index : refined[f]) ohits[hit_index * 9 + 7] = f;
    }

    ITensor::vector* tensors = new ITensor::vector;
    {
        Configuration md;
        md["name"] = "opflash";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nflash, mcol}, matrix.data(), md));
    }
    {
        Configuration md;
        md["name"] = "flash_summary";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nflash, (size_t)8}, summary.data(), md));
    }
    {
        Configuration md;
        md["name"] = "ophits";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nhit, (size_t)9}, ohits.data(), md));
    }

    Configuration md = in->metadata();
    md["producer"] = "wct-flash";
    md["nchan"] = m_nchan;
    md["offset_us"] = m_offset_us;   // per-event trigger offset for downstream Q/L
    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), md,
                                                 ITensor::shared_vector(tensors));
    log->debug("set {}: {} flashes from {} hits", in->ident(), nflash, nhit);
    return true;
}
