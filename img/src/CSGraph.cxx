#include "WireCellImg/CSGraph.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/SimpleBlob.h"

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Img::CS;



indexed_vdescs_t CS::select_ordered(const graph_t& csg,
                                    node_t::Kind kind)
{
    std::vector<vdesc_t> ret;
    for (const auto& v : vertex_range(csg)) {
        const auto& vp = csg[v];
        if (vp.kind == kind) {
            ret.push_back(v);
        }
    }
    std::sort(ret.begin(), ret.end(), [&](const vdesc_t& a, const vdesc_t& b) {
        return csg[a].ordering < csg[b].ordering;
    });
    return indexed_vdescs_t(ret);
}

graph_t CS::solve(const graph_t& csg, const SolveParams& params)
{
    graph_t csg_out;

    // Copy graph level properties
    csg_out[boost::graph_bundle] = csg[boost::graph_bundle];

    IndexedSet<vdesc_t> blob_descs_out;
    auto blob_descs = select_ordered(csg, node_t::blob);
    const size_t nblob = blob_descs.size();
    if (!nblob) {
        // std::cerr << "no blobs from graph with " << boost::num_vertices(csg) << std::endl;
        return csg_out;
    }

    IndexedSet<vdesc_t> meas_descs_out;
    auto meas_descs = select_ordered(csg, node_t::meas);
    const size_t nmeas = meas_descs.size();
    if (!nmeas) {
        // std::cerr << "no meas from graph with " << boost::num_vertices(csg) << std::endl;
        return csg_out;
    }

    double_vector_t source = double_vector_t::Zero(nblob);
    double_vector_t weight = double_vector_t::Zero(nblob);
    for (size_t bind=0; bind<nblob; ++bind) {
        const auto& blob_in = csg[blob_descs.collection[bind]];

        const auto valerr = blob_in.value;
        source(bind) = valerr.value();
        weight(bind) = valerr.uncertainty();

        auto desc_out = boost::add_vertex(blob_in, csg_out);
        blob_descs_out(desc_out);        
    }

    double_vector_t measure = double_vector_t::Zero(nmeas);
    double_matrix_t mcov = double_matrix_t::Zero(nmeas, nmeas);
    for (size_t mind=0; mind<nmeas; ++mind) {
        const auto& meas_in = csg[meas_descs.collection[mind]];
        const auto valerr = meas_in.value;
        measure(mind) = valerr.value();
        mcov(mind, mind) = valerr.uncertainty();

        auto desc_out = boost::add_vertex(meas_in, csg_out);
        meas_descs_out(desc_out);        
    }
    if (params.whiten and mcov.sum() == 0.0) {
        // std::cerr << "zero measure covariance from " << boost::num_vertices(csg) << " node graph\n";
        return csg_out;
    }

    // special case of one blob
    if (blob_descs.size() == 1) {
        auto nbdesc = blob_descs_out.collection[0];
        value_t val;
        for (size_t mind=0; mind < nmeas; ++mind) {
            auto nmdesc = meas_descs_out.collection[mind];
            boost::add_edge(nbdesc, nmdesc, csg_out);
            val += csg_out[nmdesc].value;
        }
        csg_out[nbdesc].value = val / nmeas;
        return csg_out;
    }
        
    double_matrix_t A = double_matrix_t::Zero(nmeas, nblob);

    for (auto [ei, ei_end] = boost::edges(csg); ei != ei_end; ++ei) {
        const vdesc_t tail = boost::source(*ei, csg);
        const vdesc_t head = boost::target(*ei, csg);

        const auto kind1 = csg[tail].kind;
        const auto kind2 = csg[head].kind;
        int bind{0}, mind{0};
        if (kind1 == node_t::blob and kind2 == node_t::meas) {
            bind = blob_descs.get(tail);
            mind = meas_descs.get(head);
        }        
        else if (kind2 == node_t::blob and kind1 == node_t::meas) {
            bind = blob_descs.get(head);
            mind = meas_descs.get(tail);
        }
        else {
            // someone has violated my requirements with this edge!
            continue;
        }
        A(mind, bind) = 1;

        boost::add_edge(blob_descs_out.collection[bind],
                        meas_descs_out.collection[mind],
                        csg_out);
    }

    
    double_vector_t m_vec = measure;
    double_matrix_t R_mat = A;

    if (params.whiten) {

        // std::cerr << "A:\n" << A << "\nmcov:\n" << mcov << "\nmcovinv:\n" << mcov.inverse() << std::endl;
        Eigen::LLT<double_matrix_t> llt(mcov.inverse());
        double_matrix_t U = llt.matrixL().transpose();
        // std::cerr << "U:\n" << U << std::endl;
 
        // The measure vector in a "whitened" basis
        m_vec = params.scale*U*measure;

        // The blob-measure association in "whitened" basis (becomes
        // the "reasponse" matrix in ress solving).
        R_mat = U * A;
    }

    // std::cerr << "R:\n" << R_mat << "\nm:\n" << m_vec << std::endl;
    auto solution = Ress::solve(R_mat, m_vec, params.ress,
                                source, weight);
    auto predicted = Ress::predict(R_mat, solution);

    auto& gp_out = csg_out[boost::graph_bundle];
    gp_out.chi2_base = Ress::chi2_base(m_vec, predicted);
    gp_out.chi2_l1 = Ress::chi2_l1(m_vec, solution, params.ress.lambda);

    // Update outgoing blob nodes with their solution
    for (size_t ind=0; ind<nblob; ++ind) {
        auto& bvalue = csg_out[blob_descs_out.collection[ind]];
        //bvalue.value.value(solution[ind]);
        bvalue.value = solution[ind]; // drops weight
    }
    return csg_out;
}

graph_t CS::prune(const graph_t& csg, float threshold)
{
    graph_t csg_out;

    // Copy graph level properties
    csg_out[boost::graph_bundle] = csg[boost::graph_bundle];
    
    size_t nblobs = 0;
    std::unordered_map<vdesc_t, vdesc_t> old2new;
    for (auto oldv : vertex_range(csg)) {
        const auto& node = csg[oldv];
        if (node.kind == node_t::blob) {
            if (node.value.value() < threshold) {
                continue;
            }
            ++nblobs;
        }
        old2new[oldv] = boost::add_vertex(node, csg_out);
    }
    
    if (!nblobs) {
        return csg_out;
    }

    for (auto edge : mir(boost::edges(csg))) {
        auto old_tail = boost::source(edge, csg);
        auto old_head = boost::target(edge, csg);

        auto old_tit = old2new.find(old_tail);
        if (old_tit == old2new.end()) {
            continue;
        }
        auto old_hit = old2new.find(old_head);
        if (old_hit == old2new.end()) {
            continue;
        }
        boost::add_edge(old_tit->second, old_hit->second, csg_out);
    }    
    return csg_out;
}

void CS::connected_subgraphs(const graph_t& slice_graph,
                             std::back_insert_iterator<graph_vector_t> subgraphs_out)
{
    auto islice = slice_graph[boost::graph_bundle].islice;
    if (!islice) {
        THROW(ValueError() << errmsg{"Slice graph has no slice"});
    }

    std::unordered_map<int, std::vector<vdesc_t> > groups;
    std::unordered_map<vdesc_t, int> desc2id;
    boost::connected_components(slice_graph, boost::make_assoc_property_map(desc2id));
    for (auto& [desc,id] : desc2id) {  // invert
        groups[id].push_back(desc);
    }
    for (const auto& [index, slice_vdescs] : groups) {
        graph_t sub_graph;
        sub_graph[boost::graph_bundle].islice = islice;
        sub_graph[boost::graph_bundle].index = index;
        std::unordered_map<vdesc_t, vdesc_t> slice2sub;
        // Add map from slice to sub graph vdesc's
        for (const auto& slice_vdesc : slice_vdescs) {
            vdesc_t sub_vdesc = boost::add_vertex(slice_graph[slice_vdesc], sub_graph);
            slice2sub[slice_vdesc] = sub_vdesc;
        }
        // Add edges
        for (const auto& slice_vdesc : slice_vdescs) {
            for (auto edge : mir(boost::out_edges(slice_vdesc, slice_graph))) {
                auto tail = slice2sub[boost::source(edge, slice_graph)];
                auto head = slice2sub[boost::target(edge, slice_graph)];
                boost::add_edge(tail, head, sub_graph);
            }
        }
        subgraphs_out = sub_graph;
    }
}


void CS::unpack(const cluster_graph_t& cgraph,
                std::back_insert_iterator<graph_vector_t> subgraphs,
                const value_t& meas_thresh)
{
    auto cnull = boost::graph_traits<cluster_graph_t>::null_vertex();
    auto sgnull = boost::graph_traits<graph_t>::null_vertex();

    // map slice to the slice graph
    std::unordered_map<slice_t, graph_t> s2slg;
    // map some cg vertices to slice graph vertices.  This is sparse
    // and spans the slice graphs.
    std::vector<vdesc_t> c2slg(boost::num_vertices(cgraph), sgnull);

    for (auto bvtx : vertex_range(cgraph)) {
        const auto& bnode = cgraph[bvtx];
        if (bnode.code() != 'b') {
            continue;
        }

        cluster_vertex_t svtx = cnull;
        std::vector<cluster_vertex_t> mvtxs;

        // find neighbor s-node and m-nodes
        for (auto nnvtx : mir(boost::adjacent_vertices(bvtx, cgraph))) {
            const auto& cnode = cgraph[nnvtx];
            char code = cnode.code();
            switch(code) {
                case 's': svtx=nnvtx; break;
                case 'm': mvtxs.push_back(nnvtx); break;
                default: continue;
            }
        }
        if (svtx == cnull) {
            THROW(ValueError() << errmsg{"Slice graph has no slice"});
        }
        if (mvtxs.empty()) {
            THROW(ValueError() << errmsg{"Blob has no measures"});
        }

        // get slice graph, initialize if new
        const auto& islice = std::get<slice_t>(cgraph[svtx].ptr);
        auto& slg = s2slg[islice];
        auto& slgprop = slg[boost::graph_bundle];
        if (! slgprop.islice) {
            slgprop.islice = islice;
            slgprop.index = s2slg.size() - 1;
        }

        // add new or reuse prior mnodes, applying threshold on new
        std::vector<vdesc_t> meas_descs;
        for (auto mvtx : mvtxs) {
            auto meas_desc = c2slg[mvtx];
            if (meas_desc == sgnull) {
                const auto& mnode = cgraph[mvtx];
                auto mptr = std::get<meas_t>(mnode.ptr);
                const auto msum = mptr->signal();
                if (msum.value() < meas_thresh.value() or
                    msum.uncertainty() > meas_thresh.uncertainty()) {
                    continue;
                }
                const int ordering = mnode.ident();
                node_t meas{mvtx, node_t::meas, ordering, msum};
                meas_desc = boost::add_vertex(meas, slg);
                c2slg[mvtx] = meas_desc;
            }
            meas_descs.push_back(meas_desc);
        }
        if (meas_descs.empty()) {
            continue;
        }
        
        // make b-m edges
        auto bptr = std::get<blob_t>(bnode.ptr);
        const int ordering = bnode.ident();
        node_t blob{bvtx, node_t::blob, ordering, value_t(0,1)};
        auto blob_desc = boost::add_vertex(blob, slg);

        for (auto meas_desc : meas_descs) {
            boost::add_edge(blob_desc, meas_desc, slg);
        }
    }

    // c2slg now holds per slice b-m graphs.  Each of those graphs are
    // now factored into connected subgraphs.

    for (const auto& slit : s2slg) {
        connected_subgraphs(slit.second, subgraphs);
        // subgraphs = slit.second;
    }
}

// Return new blob same as the old blob, but new value
IBlob::pointer replace_blob_value(value_t nv, cluster_vertex_t v, const cluster_graph_t& g)
{
    IBlob::pointer ob = get<IBlob::pointer>(g[v].ptr);
    return std::make_shared<Aux::SimpleBlob>(
        ob->ident(), nv.value(), nv.uncertainty(), 
        ob->shape(), ob->slice(), ob->face());
}

// The csgs provide b-m graph with vdesc into cgraph.  These {b} and
// {m} sets are subsets of those from the cgraph.  We must replace
// IBlobs on b nodes, remove b and m but leave the rest of the cgraph
// as-is.
cluster_graph_t CS::repack(const cluster_graph_t& cgin,
                           const std::vector<graph_t>& csgs)
{
    // Note, removal of vertices from cluster_graph_t is dreadfully
    // slow due to using vecS to hold vertices.  Thus, we use a
    // somewhat verbose construction instead of copy+remove.

    // Map old blob vertex to new blob value.
    std::unordered_map<cluster_vertex_t, value_t> live_bs;
    // Just remember which m's survive
    std::unordered_set<cluster_vertex_t> live_ms;

    // Pass over result graphs to index the live-{b,m}'s by their
    // vertex descriptors of cgin.
    for (const graph_t& csg : csgs) {
        for (auto v : vertex_range(csg)) {
            const auto& node = csg[v];
            if (node.kind == node_t::meas) {
                live_ms.insert(node.orig_desc);
                continue;
            }
            if (node.kind == node_t::blob) {
                live_bs[node.orig_desc] = node.value;
                // std::cerr << "blob " << node.orig_desc
                //           << " q=" << node.value << "\n";
                continue;
            }
        }
    }

    // Pass over input graph and copy any non-{b,m} and any
    // live-{b,m}.
    cluster_graph_t cgtmp;
    std::unordered_map<cluster_vertex_t, cluster_vertex_t> old2new;
    for (auto oldv : vertex_range(cgin)) {
        const auto& node = cgin[oldv];
        const char code = node.code();
        if (code == 'b') {
            auto it = live_bs.find(oldv);
            if (it == live_bs.end()) {
                continue;
            }
            auto newb = replace_blob_value(it->second, oldv, cgin);
            auto newv = boost::add_vertex(newb, cgtmp);
            old2new[oldv] = newv;
            continue;
        }

        if (code == 'm') {
            if (live_ms.count(oldv)) {
                auto newv = boost::add_vertex(node, cgtmp);
                old2new[oldv] = newv;
            }
            continue;
        }
        // o.w. its a non-{b,m}
        auto newv = boost::add_vertex(node, cgtmp);
        old2new[oldv] = newv;

        // debug
        // if (code == 's') {
        //     auto islice = get<ISlice::pointer>(node.ptr);
        //     ISlice::value_t tot;
        //     for (const auto& [ich, val] : islice->activity()) {
        //         tot += val;
        //     }
        //     std::cerr << "slice " << islice->ident() << " " << tot << std::endl;
        // }
    }

    // Pass over input graph edges and transfer any where we know both
    // ends.
    for (auto edge : mir(boost::edges(cgin))) {
        auto old_tail = boost::source(edge, cgin);
        auto old_head = boost::target(edge, cgin);

        auto old_tit = old2new.find(old_tail);
        if (old_tit == old2new.end()) {
            continue;
        }
        auto old_hit = old2new.find(old_head);
        if (old_hit == old2new.end()) {
            continue;
        }
        boost::add_edge(old_tit->second, old_hit->second, cgtmp);
    }

    // At this point the graph may have neighborless b's and m's.  We
    // could do another pass on cgtmp to construct a final out, but
    // for now leave these orphans as-is.
    return cgtmp;
}

