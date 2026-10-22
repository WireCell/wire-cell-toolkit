// Duplicated from img/src/BlobDepoFill.cxx (2026-10-21, toolkit ab262c28) with
// the per-track accumulation added.  Keep the integration itself in step with
// the original: the parity doctest and the gate on real events assert that
// filled_cluster() is bit-identical to BlobDepoFill.

#include "WireCellTruth/DepoBlobAssoc.h"

#include "WireCellAux/SliceTools.h"
#include "WireCellAux/DepoTools.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleCluster.h"

#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"

#include <boost/graph/copy.hpp>
#include <cmath>
#include <limits>
#include <unordered_set>

using namespace WireCell;
using WireCell::Aux::sensitive;
using WireCell::GraphTools::mir;

// A directed graph holding s-d-w nodes with weighted edges giving the
// integral of the Gausian section.  Always add edges in order as:
// add_edge(s,d) and add_edge(d,w)
namespace idweight {
    struct node_t { size_t idx; char code; };
    struct edge_t { double weight; };
    using graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS,
                                          node_t, edge_t>;
    using vdesc_t = boost::graph_traits<graph_t>::vertex_descriptor;
    using edesc_t = boost::graph_traits<graph_t>::edge_descriptor;
}

static double smeared(double sigma, double extra)
{
    // Only touch the number when the knob is on so that the off path
    // produces the very same doubles as BlobDepoFill.
    if (extra > 0.0) {
        return std::hypot(sigma, extra);
    }
    return sigma;
}

// Produce an s-d-w graph with nodes packed in [{s},{d},{w}] order
static
idweight::graph_t
slice_and_dice_depos(const Binning& sbins,       // run over s nodes
                     const IDepo::vector& depos, // run over d nodes
                     const std::vector<double>& tran_sigma, // per depo, smeared
                     const Pimpos& pimpos,       // run over w nodes
                     const Truth::DepoBlobAssocParams& params)
{
    idweight::graph_t gr;

    // Add slice nodes
    const size_t nsbins = sbins.nbins();
    for (size_t idx = 0; idx < nsbins; ++idx) {
        const auto gidx [[maybe_unused]] = boost::add_vertex({idx, 's'}, gr);
        assert(gidx == idx);
    }

    // Add depo nodes
    const size_t ndepos = depos.size();
    for (size_t idx=0; idx<ndepos; ++idx) {
        const auto gidx [[maybe_unused]] = boost::add_vertex({idx, 'd'}, gr);
        assert(gidx == idx + nsbins);
    }

    // Add wire nodes
    const auto pbins = pimpos.region_binning();
    const size_t nwires = pbins.nbins();
    for (size_t idx = 0; idx < nwires; ++idx) {
        const auto gidx [[maybe_unused]] = boost::add_vertex({idx, 'w'}, gr);
        assert(gidx == idx + nsbins + ndepos);
    }

    // Make s-d and d-w edges holding Gaussian integration values.
    for (size_t depo_idx=0; depo_idx<ndepos; ++depo_idx) {
        const auto& idepo = depos[depo_idx];
        const auto ddesc = depo_idx + nsbins;

        // first we slice
        {
            const double tmean = idepo->time() + params.time_offset;
            const double tsigma = smeared(idepo->extent_long()/params.speed, params.extra_smear_long);
            const double tmin = tmean - tsigma*params.nsigma;
            const double tmax = tmean + tsigma*params.nsigma;
            const Binning tbins = subset(sbins, tmin, tmax);
            const auto ntbins = tbins.nbins();
            if (ntbins) {
                std::vector<double> weights(ntbins);
                gaussian(weights.begin(), tbins, tmean, tsigma);

                for (int it=0; it<ntbins; ++it) {
                    const idweight::vdesc_t sdesc = sbins.bin(tbins.center(it));
                    boost::add_edge(sdesc, ddesc, {weights[it]}, gr);
                }
            }
        }

        // then we dice
        {
            // (x,y,z) -> (normal, wire, pitch);
            const auto wpos = pimpos.transform(idepo->pos());
            const double wmean = wpos[2]; // pitch
            const double wsigma = tran_sigma[depo_idx];

            auto wbins = subset(pbins, wmean-wsigma*params.nsigma, wmean+wsigma*params.nsigma);
            const auto nwbins = wbins.nbins();
            if (nwbins) {
                std::vector<double> weights(wbins.nbins(), 0);
                gaussian(weights.begin(), wbins, wmean, wsigma);

                // make edge from each d to w
                int wip = pbins.bin(wbins.center(0));
                for (int wind = 0; wind < wbins.nbins(); ++wind) {
                    const idweight::vdesc_t wdesc = wip + wind + nsbins + ndepos;
                    boost::add_edge(ddesc, wdesc, {weights[wind]}, gr);
                }
            }
        }
    }

    return gr;
}

ICluster::pointer Truth::filled_cluster(const ICluster::pointer& icluster,
                                        const DepoBlobAssocResult& res)
{
    cluster_graph_t gr;
    boost::copy_graph(icluster->graph(), gr);

    // Collect the info from just b-type nodes
    for (const auto& vdesc : mir(boost::vertices(gr))) {
        const char code = gr[vdesc].code();
        if (code == 'b') {
            auto bold = std::get<IBlob::pointer>(gr[vdesc].ptr);
            double val = 0.0;
            auto it = res.blobs.find(vdesc);
            if (it != res.blobs.end()) {
                val = it->second.q_total;
            }
            // make a new IBlob with new charge and set it as payload
            auto bnew = std::make_shared<Aux::SimpleBlob>(
                bold->ident(), val, 0.0,
                bold->shape(), bold->slice(), bold->face());
            gr[vdesc].ptr = bnew;
        }
    }
    return std::make_shared<Aux::SimpleCluster>(std::move(gr), icluster->ident());
}

Truth::DepoBlobAssocResult
Truth::depo_blob_assoc(const ICluster::pointer& icluster,
                       const IDepo::vector& alldepos,
                       const DepoBlobAssocParams& params,
                       Log::logptr_t log)
{
    if (params.speed <= 0.0) {
        THROW(ValueError() << errmsg{"DepoBlobAssoc: no speed provided"});
    }
    if (!log) {
        log = Log::logger("truth");
    }

    DepoBlobAssocResult res;
    const auto& ingr = icluster->graph();

    // Form sets of blobs with common anode face.
    std::unordered_map<IAnodeFace::pointer, std::vector<cluster_vertex_t>> face2blobs;

    // Collect the info from just b-type nodes
    for (const auto& vdesc : mir(boost::vertices(ingr))) {
        const char code = ingr[vdesc].code();
        if (code == 'b') {
            auto iblob = std::get<IBlob::pointer>(ingr[vdesc].ptr);
            face2blobs[iblob->face()].push_back(vdesc);
            res.blobs[vdesc] = {};
        }
    }
    res.nblobs = res.blobs.size();
    if (res.blobs.empty()) {
        return res;
    }

    // Now, Process each face separately.
    for (const auto& [ianodeface, bdescvector] : face2blobs) {

        auto iwireplane = ianodeface->planes()[params.pindex];
        const Pimpos* pimpos = iwireplane->pimpos();
        const auto& coords = ianodeface->raygrid();

        // Collect the slices spanned by this face's blobs and map
        // each them to the graph vertex descriptors of the blobs in
        // the slice.
        std::unordered_set<ISlice::pointer> slices;
        std::unordered_map<ISlice::pointer, std::vector<cluster_vertex_t>> slice2blobs;
        for (const auto& bdesc : bdescvector) {
            auto iblob = std::get<IBlob::pointer>(ingr[bdesc].ptr);
            auto islice = iblob->slice();
            slices.insert(islice);
            slice2blobs[islice].push_back(bdesc);
        }

        // Make a Binning that spans the slices' start and span times.
        const auto sbins = Aux::binning(slices.begin(), slices.end());

        // Get depos in the sensitive area of the current face and
        // pre-calculate their central position in the primary wire
        // plane "pimpos" coordinate system, their track key and the
        // (smeared) transverse sigma.
        IDepo::vector depos;
        std::vector<double> dwcenter, tran_sigma, key_energy;
        std::vector<int> key_id;
        auto mydepos = sensitive(alldepos, ianodeface);
        const int dirx = ianodeface->dirx();
        const double behind_tolerance = 0.01*units::mm;
        for (const auto& maybe : mydepos) {
            const auto rpos = pimpos->relative(maybe->pos());
            if (dirx * rpos[0] < -behind_tolerance) {
                continue;       // depo is behind the face.
            }
            depos.push_back(maybe);
            res.q_depos += std::abs(maybe->charge());
            dwcenter.push_back(pimpos->transform(maybe->pos())[1]);
            tran_sigma.push_back(smeared(maybe->extent_tran(), params.extra_smear_tran));
            const IDepo* src = maybe.get();
            if (params.id_from_prior && maybe->prior()) {
                src = maybe->prior().get();
            }
            int tid = src->id();
            double energy = src->energy();
            if (params.energy_map && tid >= 0 && (size_t)tid < params.energy_map->size()) {
                energy = (*params.energy_map)[tid];
            }
            if (params.id_map && tid >= 0 && (size_t)tid < params.id_map->size()) {
                tid = (*params.id_map)[tid];
            }
            key_id.push_back(tid);
            key_energy.push_back(energy);
        }
        res.ndepos += depos.size();
        for (const auto& d : depos) {
            const double t = d->time() + params.time_offset;
            if (t >= sbins.min() && t < sbins.max()) {
                ++res.ndepos_window;
                res.q_depos_window += std::abs(d->charge());
            }
        }
        log->debug("depo_blob_assoc: face={} nblobs={} ndepos={} in window {}",
                   ianodeface->ident(), bdescvector.size(), depos.size(), res.ndepos_window);

        auto gr = slice_and_dice_depos(sbins, depos, tran_sigma, *pimpos, params);

        // Per-slice blob-by-wire index, as BlobDepoFill.
        struct blobwire_t {
            DepoBlobAssocResult::BlobEntry* entry;
            double wlo, whi;    // blob extent along the primary wire, from the other planes
        };
        const int nfake = 2;
        const int nplanes = coords.nlayers() - nfake;

        // Each slice
        for (const auto& [islice,bdescs] : slice2blobs) {
            const auto sdesc = sbins.bin(islice->start() + 0.5*islice->span());

            // Index the slice's blobs by primary-plane wire, in bdescs order.
            std::unordered_map<int, std::vector<blobwire_t>> bywire;
            for (const auto& bdesc : bdescs) {
                const IBlob* iblob = std::get<IBlob::pointer>(ingr[bdesc].ptr).get();
                const auto& strips = iblob->shape().strips();
                const auto& wiprange = strips[2 + params.pindex].bounds;
                auto* entry = &res.blobs[bdesc];
                for (int wip = wiprange.first; wip < wiprange.second; ++wip) {
                    RayGrid::coordinate_t rgc{2+params.pindex, wip};
                    double wlo = -std::numeric_limits<double>::infinity();
                    double whi = std::numeric_limits<double>::infinity();
                    for (int irel=1; irel<nplanes; ++irel) {
                        int other_pindex = (params.pindex + irel)%nplanes;
                        const auto owipr = strips[nfake + other_pindex].bounds;

                        RayGrid::coordinate_t orgc1{nfake+other_pindex, owipr.first};
                        RayGrid::coordinate_t orgc2{nfake+other_pindex, owipr.second};

                        auto p1 = coords.ray_crossing(rgc, orgc1);
                        auto p2 = coords.ray_crossing(rgc, orgc2);

                        auto w1 = pimpos->transform(p1)[1];
                        auto w2 = pimpos->transform(p2)[1];
                        if (w2 < w1) std::swap(w1,w2);

                        wlo = std::max(wlo, w1);
                        whi = std::min(whi, w2);
                    }
                    bywire[wip].push_back({entry, wlo, whi});
                }
            }

            // Each depo in the slice
            for (auto sdedge : mir(boost::out_edges(sdesc, gr))) {
                const double sd_weight = gr[sdedge].weight;
                const auto ddesc = boost::target(sdedge, gr);
                const auto didx = gr[ddesc].idx;
                const auto& idepo = depos[didx];
                const double depo_center = dwcenter[didx];
                const double depo_sigma = tran_sigma[didx];
                const double depo_min = depo_center - depo_sigma*params.nsigma;
                const double depo_max = depo_center + depo_sigma*params.nsigma;
                const double depo_q = std::abs(idepo->charge());
                const int tid = key_id[didx];
                const double energy = key_energy[didx];

                // Each wire hit by the depo
                for (auto dwedge : mir(boost::out_edges(ddesc, gr))) {
                    const double dw_weight = gr[dwedge].weight;
                    const auto wdesc = boost::target(dwedge, gr);
                    const int wip = gr[wdesc].idx;

                    const auto bwit = bywire.find(wip);
                    if (bwit == bywire.end()) {
                        continue;   // no blob of the slice on this wire
                    }
                    for (const auto& bw : bwit->second) {
                        const double w1 = std::max(bw.wlo, depo_min);
                        const double w2 = std::min(bw.whi, depo_max);
                        if (w1 >= w2) {
                            continue; // no overlap (the wcfm doc 05 defect)
                        }

                        const double w_weight = gbounds(w1,w2,depo_center, depo_sigma);
                        const double dq = sd_weight * dw_weight * w_weight * depo_q;
                        bw.entry->q_total += dq;
                        auto& bt = bw.entry->by_track[tid];
                        bt.first += dq;
                        if (depo_q > 0.0) {
                            bt.second += energy * (dq / depo_q);
                        }
                        if (params.keep_depos) {
                            auto& ds = bw.entry->depos;
                            if (!ds.empty() && ds.back().depo == idepo.get()) {
                                ds.back().dq += dq;   // the same depo on the next wire
                            }
                            else {
                                ds.push_back({idepo.get(), tid, dq});
                            }
                        }
                    } // over blobs on the wire
                } // over wires
            } // over s-d edges
        } // over slices
    } // over faces

    return res;
}
