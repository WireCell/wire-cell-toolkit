#include "WireCellImg/BlobDepoFill.h"

#include "WireCellAux/SliceTools.h"
#include "WireCellAux/DepoTools.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/ClusterHelpers.h"

#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/NamedFactory.h"


WIRECELL_FACTORY(BlobDepoFill, WireCell::Img::BlobDepoFill,
                 WireCell::INamed,
                 WireCell::IBlobDepoFill, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::sensitive;
using WireCell::GraphTools::mir;

Img::BlobDepoFill::BlobDepoFill()
    : Aux::Logger("BlobDepoFill", "img")
{
}
Img::BlobDepoFill::~BlobDepoFill()
{
}

WireCell::Configuration Img::BlobDepoFill::default_configuration() const
{
    Configuration cfg;

    // When drifted depos are sent through DepoTransform to form the
    // input to blobs, they gain a time offset due to the "start_time"
    // config parameter.  This provided time offset is added to depo
    // times prior to associating with blobs.
    cfg["time_offset"] = m_toffset;

    // The number of sigma of depo Gaussian to consider.
    cfg["nsigma"] = m_nsigma;
    
    // The nominal drift speed used to convert the longitudinal extent
    // of depos from units [distance] to units of [time].
    cfg["speed"] = m_speed;

    // The plane index to use to use as the "primary" plane.  The
    // depos will be "diced" onto the wires of this plane.
    cfg["pindex"] = m_pindex;

    return cfg;
}


void Img::BlobDepoFill::configure(const WireCell::Configuration& cfg)
{
    m_toffset = get(cfg, "time_offset", m_toffset);
    m_nsigma = get(cfg, "nsigma", m_nsigma);
    m_speed = get(cfg, "speed", m_speed);
    m_pindex = get(cfg, "pindex", m_pindex);

    if (m_speed <= 0.0) {
        THROW(ValueError() << errmsg{"no speed provided"});
    }
}

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

// Produce an s-d-w graph with nodes packed in [{s},{d},{w}] order
static
idweight::graph_t
slice_and_dice_depos(Log::logptr_t& log,
                     const Binning& sbins,       // run over s nodes
                     const IDepo::vector& depos, // run over d nodes
                     const Pimpos& pimpos,       // run over w nodes
                     double speed, // to convert longitudinal extent to time
                     double nsigma, double time_offset)
{
    idweight::graph_t gr;
    
    // Add slice nodes 
    const size_t nsbins = sbins.nbins();
    for (size_t idx = 0; idx < nsbins; ++idx) {
        const auto gidx = boost::add_vertex({idx, 's'}, gr);
        assert(gidx == idx);
    }

    // Add depo nodes
    const size_t ndepos = depos.size();
    for (size_t idx=0; idx<ndepos; ++idx) {
        const auto gidx = boost::add_vertex({idx, 'd'}, gr);
        assert(gidx == idx + nsbins);
    }
    
    // Add wire nodes
    const auto pbins = pimpos.region_binning();
    const size_t nwires = pbins.nbins();
    for (size_t idx = 0; idx < nwires; ++idx) {
        const auto gidx = boost::add_vertex({idx, 'w'}, gr);
        assert(gidx == idx + nsbins + ndepos);
    }

    // Make s-d and d-w edges holding Gaussian integration values.
    for (size_t depo_idx=0; depo_idx<ndepos; ++depo_idx) {
        const auto& idepo = depos[depo_idx];
        const auto ddesc = depo_idx + nsbins;

        // first we slice
        {
            const double tmean = idepo->time() + time_offset;
            const double tsigma = idepo->extent_long()/speed;
            const double tmin = tmean - tsigma*nsigma;
            const double tmax = tmean + tsigma*nsigma;
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
            const double wsigma = idepo->extent_tran();

            auto wbins = subset(pbins, wmean-wsigma*nsigma, wmean+wsigma*nsigma);
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

static 
ICluster::pointer make_cluster(const ICluster::pointer& icluster,
                        std::unordered_map<cluster_vertex_t, std::pair<double, double>>& blob_value)
{
    cluster_graph_t gr;
    boost::copy_graph(icluster->graph(), gr);

    // Collect the info from just b-type nodes
    for (const auto& vdesc : mir(boost::vertices(gr))) {
        const char code = gr[vdesc].code();
        if (code == 'b') {
            auto bold = std::get<IBlob::pointer>(gr[vdesc].ptr);
            // make a new IBlob with new charge and set it as payload
            auto [val, unc] = blob_value[vdesc];
            auto bnew = std::make_shared<Aux::SimpleBlob>(
                bold->ident(), val, unc,
                bold->shape(), bold->slice(), bold->face());
            gr[vdesc].ptr = bnew;
        }
    }
    return std::make_shared<Aux::SimpleCluster>(gr, icluster->ident());
}


bool Img::BlobDepoFill::operator()(const input_tuple_type& intup,
                                   output_pointer& out)
{
    out = nullptr;
    auto icluster = std::get<0>(intup);
    auto ideposet = std::get<1>(intup);
    if (!icluster or !ideposet) {
        // If one is EOS, both should be EOS
        if (icluster) {
            const auto blobs = oftype<IBlob::pointer>(icluster->graph());
            log->warn("join node got non-EOS ICluster {} with {} blobs on EOS",
                      icluster->ident(), blobs.size());
        }
        if (ideposet) {
            log->warn("join node got non-EOS IDepoSet {} with {} depos on EOS",
                      ideposet->ident(),
                      ideposet->depos()->size());
        }
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    log->debug("call={} fill cluster {} with deposet {}",
               m_count, icluster->ident(), ideposet->ident());

    const auto& ingr = icluster->graph();

    // We will collect the new blob charge by its vertex descriptor.
    std::unordered_map<cluster_vertex_t, std::pair<double, double>> blob_value;

    // Form sets of blobs with common anode face.
    std::unordered_map<IAnodeFace::pointer, std::vector<cluster_vertex_t>> face2blobs;

    // Collect the info from just b-type nodes
    for (const auto& vdesc : mir(boost::vertices(ingr))) {
        const char code = ingr[vdesc].code();
        if (code == 'b') {
            auto iblob = std::get<IBlob::pointer>(ingr[vdesc].ptr);
            face2blobs[iblob->face()].push_back(vdesc);
            blob_value[vdesc] = {0.0, 0.0};
        }
    }

    if (blob_value.empty()) {
        out = make_cluster(icluster, blob_value);
        log->debug("call={} no input blobs", m_count);
        ++m_count;
        return true;
    }

    // just for logging.
    size_t nblobs=0;

    // Now, Process each face separately.
    for (const auto& [ianodeface, bdescvector] : face2blobs) {

        auto iwireplane = ianodeface->planes()[m_pindex];
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

        // Make a Binning htat spans the slices' start and span times.
        // Because we take care to fill the s-d-w graph below in a
        // particular order, the sbins indices are also directly the
        // graphs's s-node descriptors.
        const auto sbins = Aux::binning(slices.begin(), slices.end());

        // Get depos in the sensitive area of the current face and
        // pre-calculate their central position in the primary wire
        // plane "pimpos" coordinate system.
        IDepo::vector depos;
        std::vector<double> dwcenter;
        auto mydepos = sensitive(*(ideposet->depos()), ianodeface);
        log->debug("call={} face={} nblobs={} ndepos={}",
                   m_count, bdescvector.size(), mydepos.size());
        for (const auto& maybe : mydepos) {
            const auto rpos = pimpos->relative(maybe->pos());
            if (rpos[0] < 0) {
                continue;       // depo is behind the face.
            }
            depos.push_back(maybe);
            dwcenter.push_back(pimpos->transform(maybe->pos())[1]);
        }

        // Make the s-d-w graph holding Gaussian integrals over slice
        // and primary wire directions.  Note, nodes are indexed in
        // order to keep each type contiguous: [{s},{d},{w}].  {s} is
        // ordered as per slice Binning, {d} ordered per deposet, {w}
        // ordered per pimpos wires.
        auto gr = slice_and_dice_depos(log, sbins, depos, *pimpos,
                                       m_speed, m_nsigma, m_toffset);

        // Each slice
        for (const auto& [islice,bdescs] : slice2blobs) {
            const auto sdesc = sbins.bin(islice->start() + 0.5*islice->span());
            nblobs += bdescs.size();

            // Each depo in the slice
            for (auto sdedge : mir(boost::out_edges(sdesc, gr))) {
                const double sd_weight = gr[sdedge].weight;
                const auto ddesc = boost::target(sdedge, gr);
                const auto didx = gr[ddesc].idx;
                const auto& idepo = depos[didx];
                const double depo_center = dwcenter[didx];
                const double depo_min = depo_center - idepo->extent_tran()*m_nsigma;
                const double depo_max = depo_center + idepo->extent_tran()*m_nsigma;

                // Each wire hit by the depo
                for (auto dwedge : mir(boost::out_edges(ddesc, gr))) {
                    const double dw_weight = gr[dwedge].weight;
                    const auto wdesc = boost::target(dwedge, gr);
                    const int wip = gr[wdesc].idx;

                    // Each blob in slice
                    for (const auto& bdesc : bdescs) {
                        auto iblob = std::get<IBlob::pointer>(ingr[bdesc].ptr);
        
                        // blob bounds
                        const auto& strips = iblob->shape().strips();

                        // wire in plane range for blob
                        const auto& wiprange = strips[2 + m_pindex].bounds;

                        if (wip < wiprange.first or wip >= wiprange.second) {
                            continue; // wire not in blob.
                        }

                        RayGrid::coordinate_t rgc{2+m_pindex, wip};

                        std::vector<double> lo, hi;
                        // fixme: assumes first 2 layers are "fake
                        // planes" given the overall sensitive bounds.
                        // Should always hold, but someone someday may
                        // get "creative".
                        const int nfake = 2;
                        const int nplanes = coords.nlayers() - nfake;
                        for (int irel=1; irel<nplanes; ++irel) {
                            int other_pindex = (m_pindex + irel)%nplanes;
                            const auto owipr = strips[nfake + other_pindex].bounds;

                            RayGrid::coordinate_t orgc1{nfake+other_pindex, owipr.first};
                            RayGrid::coordinate_t orgc2{nfake+other_pindex, owipr.second};

                            auto p1 = coords.ray_crossing(rgc, orgc1);
                            auto p2 = coords.ray_crossing(rgc, orgc2);

                            auto w1 = pimpos->transform(p1)[1];
                            auto w2 = pimpos->transform(p2)[1];
                            if (w2 < w1) std::swap(w1,w2);
                            
                            lo.push_back(w1);
                            hi.push_back(w2);
                        }
                        const double wlo = *std::max_element(lo.begin(), lo.end());
                        const double whi = *std::min_element(hi.begin(), hi.end());

                        const double w1 = std::max(wlo, depo_min);
                        const double w2 = std::min(whi, depo_max);

                        const double w_weight = gbounds(w1,w2,depo_center, idepo->extent_tran());
                        const double dq = sd_weight * dw_weight * w_weight * std::abs(idepo->charge());
                        auto& bv = blob_value[bdesc];
                        bv.first += dq;
                        bv.second += 0; // how to estimate?
                    } // over blobs
                } // over wires
            } // over s-d edges
        } // over slices
    } // over faces

    out = make_cluster(icluster, blob_value);
    log->debug("call={} nblobs={}", m_count, nblobs);
    ++m_count;
    return true;
}


