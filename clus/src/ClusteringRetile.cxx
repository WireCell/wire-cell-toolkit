// Developers: see important comments in header file.
//
// A "FIXME" really must be fixed if you expect anything to work.  Lower case
// "fixme" are "normal" fixmes that everyone ignores.


#include "WireCellClus/ClusteringRetile.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/RayHelpers.h"

using namespace WireCell;

// Segregate this weird choice for namespace.
namespace WCC = WireCell::PointCloud::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;


WCC::ClusteringRetile::ClusteringRetile(const WireCell::Configuration& cfg)
{
    auto sampler = get<std::string>(cfg, "sampler","");
    if (sampler.empty()) {
        raise<ValueError>("ClusteringRetile requires an IBlobSampler type/name in 'sampler' parameter");
    }
    m_sampler = Factory::find_tn<IBlobSampler>(sampler);

    auto anode_tn = get<std::string>(cfg, "anode","");
    if (anode_tn.empty()) {
        raise<ValueError>("ClusteringRetile requires an IAnodePlane type/name in 'anode' parameter");
    }
    int face_index = get(cfg, "face", 0);

    auto anode = Factory::find_tn<IAnodePlane>(anode_tn);
    m_face = anode->faces()[face_index];
    if (!m_face) {
        raise<ValueError>("ClusteringRetile got null IAnodeFace at index=%d from %s", face_index, anode_tn);
    }

    const auto& coords = m_face->raygrid();
    if (coords.nlayers() != 5) {
        raise<ValueError>("unexpected number of ray grid layers: %d", coords.nlayers());
    }

    // Add time cut configuration
    m_cut_time_low = get(cfg, "cut_time_low", -1e9);
    m_cut_time_high = get(cfg, "cut_time_high", 1e9);
}


// Step 1. Build activities from blobs in a cluster.
WRG::activities_t WCC::ClusteringRetile::get_activity(const Cluster& cluster) const
{
    const int nlayers = 2+3;
    std::vector<WRG::measure_t> measures(nlayers);

    // checkme: this assumes "iend" is the usual one-past-last aka [ibeg,iend)
    // forms a half-open range.  I'm not sure if PointTreeBuilding is following
    // this or not.

    int (WCC::Blob::*wmin[])(void) const = {
        &WCC::Blob::u_wire_index_min,
        &WCC::Blob::v_wire_index_min,
        &WCC::Blob::w_wire_index_min
    };

    int (WCC::Blob::*wmax[])(void) const = {
        &WCC::Blob::u_wire_index_max,
        &WCC::Blob::v_wire_index_max,
        &WCC::Blob::w_wire_index_max
    };
        
    const double hit=1.0;       // actual charge value does not matter to tiling.

    // what to do the first two views???
    measures[0].push_back(1);
    measures[1].push_back(1);

    // the three views ...
    for (int index=0; index<3; ++index) {
        const int layer = index + 2;
        WRG::measure_t& m = measures[layer];

        // Make each "wire" in each blob's bounds of this plane "hit".
        for (const auto* fblob : cluster.children()) {
            int ibeg = (fblob->*wmin[index])();
            int iend = (fblob->*wmax[index])();
            m.reserve(iend);
            while (ibeg < iend) {
                m[ibeg++] = hit;
            }
            //std::cout << ibeg << " " << iend << " " << index << " " << hit << std::endl;
        }
    }

    return RayGrid::make_activities(m_face->raygrid(), measures);
}


// Step 2. Modify activity to suit.
WRG::activities_t WCC::ClusteringRetile::hack_activity(const WRG::activities_t& activity) const
{
    WRG::activities_t ret;


    // FIXME: Xin, delete this line and add your "hacks".  Note, your
    // "hackgorithm" may require more arguments to come in to this method that I
    // write here.  Since I don't know what you need, I only give the minimal.
    // More can be added as needed.
    ret = activity;


    return ret;
}



// Step 3. Form IBlobs from activities.
std::vector<IBlob::pointer> WCC::ClusteringRetile::make_iblobs(const WRG::activities_t& activity) const
{
    std::vector<IBlob::pointer> ret;

    const auto& coords = m_face->raygrid();

    // Do the actual tiling.
    auto bshapes = WRG::make_blobs(coords, activity);

    std::cout << bshapes.size() << " " << activity.size() << std::endl;

    // Convert RayGrid blob shapes into IBlobs 
    const float blob_value = 0.0;  // tiling doesn't consider particular charge
    const float blob_error = 0.0;  // tiling doesn't consider particular charge
    int blob_ident=0;
    for (const auto& bshape : bshapes) {
        ISlice::pointer slice = nullptr; // fixme: okay?
        IBlob::pointer iblob = std::make_shared<Aux::SimpleBlob>(blob_ident++, blob_value,
                                                                 blob_error, bshape, slice, m_face);
        std::cout << "Test: " << iblob << std::endl;
        // FIXME: (maybe okay?) GridTiling produces an IBlobSet here which holds
        // ISlice info.  Are we losing anything important not including that
        // info?
        ret.push_back(iblob);
    }
    return ret;
}



void WCC::ClusteringRetile::operator()(WCC::Grouping& original, WCC::Grouping& shadow, cluster_set_t&) const
{
    // FIXME: With #377 fixed, we would make the shadow grouping here from
    // scratch and add it to the input map by name, eg "shadow".  Instead, we
    // smash whatever is in the 2nd grouping to fill with "shadow" clusters.  In
    // other cluster functions this second grouping is interpreted as holding
    // "dead" clusters.
    shadow.local_pcs().clear();

    for (auto* orig_cluster : original.children()) {

        // find the flash time:
        auto flash = orig_cluster->get_flash();
        int nblobs = orig_cluster->kd_blobs().size();
        
        // Apply time cut
        if (flash) {
            double flash_time = flash.time();
            if (flash_time >= m_cut_time_low && flash_time <= m_cut_time_high) {
                // std::cout << "Tests: " << nblobs << " at time " << flash_time/units::us << " " << m_cut_time_low/units::us << " " << m_cut_time_high/units::us << "\n";
                
                // get the span of indices
                auto cc = orig_cluster->get_pcarray("isolated", "perblob");
                // convert span to vector
                std::vector<int> cc_vec(cc.begin(), cc.end());
                // for (const auto& val : cc_vec) {
                //     std::cout << val << " ";
                // }
                // std::cout << std::endl;

                // use the vector for separate()
                // origi_cluster still have the original main cluster ... 
                auto splits = original.separate(orig_cluster, cc_vec);
                
                std::map<int, Cluster*> map_id_cluster = splits;
                map_id_cluster[-1] = orig_cluster;

                Cluster *shadow_orig_cluster;
                std::map<int, Cluster*> shadow_splits;

                for (auto& [id, cluster] : map_id_cluster) {
                    // make a shadow cluster
                    auto& shad_cluster = shadow.make_child();
                    shad_cluster.set_ident(cluster->ident());
                    //std::cout << cluster->ident() << std::endl;

                    if (id==-1) shadow_orig_cluster = &shad_cluster;
                    else shadow_splits[id] = &shad_cluster;

                    
                    // find the highest and lowest points
                    std::pair<geo_point_t, geo_point_t> pair_points = cluster->get_highest_lowest_points();
                    //std::cout << pair_points.first << " " << pair_points.second << std::endl;
                    int high_idx = cluster->get_closest_point_index(pair_points.first);
                    int low_idx = cluster->get_closest_point_index(pair_points.second);
                    cluster->dijkstra_shortest_paths(high_idx, false);
                    cluster->cal_shortest_path(low_idx);

                    // check path ... 
                    auto path_wcps = cluster->get_path_wcps();
                    for (const auto& wcp : path_wcps) {
                        auto point = cluster->point3d(wcp);
                        std::cout << point << std::endl;
                    }

                    // Step 1.
                    auto orig_activity = get_activity(*cluster);

                    // Step 2.
                    auto shad_activity = hack_activity(orig_activity); // may need more args

                    // Step 3.  Must make IBlobs for this is what the sampler takes.
                    auto shad_iblobs = make_iblobs(shad_activity); // may need more args

                    // Steps 4-6.
                    auto niblobs = shad_iblobs.size();
                    // Forgive me (and small-f fixme), but this is now the 3rd generation of
                    // copy-paste.  Gen 2 is in UbooneClusterSource.  OG is in
                    // PointTreeBuilding.  The reason for the copy-pastes is insufficient
                    // factoring of the de-factor standard sampling code in PointTreeBuilding.
                    // Over time, it is almost guaranteed these copy-pastes become out-of-sync.  
                    for (size_t bind=0; bind<niblobs; ++bind) {
                        if (!m_sampler) {
                            shad_cluster.make_child();
                            continue;
                        }
                        const IBlob::pointer iblob = shad_iblobs[bind];

                        // Sample the iblob, make a new blob node.
                        PointCloud::Tree::named_pointclouds_t pcs;

                        auto [pc3d, aux] = m_sampler->sample_blob(iblob, bind);
                        
                        pcs.emplace("3d", pc3d);
                        /// These seem unused and bring in yet more copy-paste code
                        // pcs.emplace("2dp0", make2dds(pc3d, angle_u));
                        // pcs.emplace("2dp1", make2dds(pc3d, angle_v));
                        // pcs.emplace("2dp2", make2dds(pc3d, angle_w));
                        const Point center = WireCell::Aux::calc_blob_center(pcs["3d"]);
                        auto scalar_ds = WireCell::Aux::make_scalar_dataset(iblob, center, pcs["3d"].get("x")->size_major(), 500*units::ns);
                        int max_wire_interval = aux.get("max_wire_interval")->elements<int>()[0];
                        int min_wire_interval = aux.get("min_wire_interval")->elements<int>()[0];
                        int max_wire_type = aux.get("max_wire_type")->elements<int>()[0];
                        int min_wire_type = aux.get("min_wire_type")->elements<int>()[0];
                        scalar_ds.add("max_wire_interval", Array({(int)max_wire_interval}));
                        scalar_ds.add("min_wire_interval", Array({(int)min_wire_interval}));
                        scalar_ds.add("max_wire_type", Array({(int)max_wire_type}));
                        scalar_ds.add("min_wire_type", Array({(int)min_wire_type}));
                        pcs.emplace("scalar", std::move(scalar_ds));

                        shad_cluster.node()->insert(Tree::Points(std::move(pcs)));
                        std::cout << "test " << std::endl;
                    }

                    std::cout << m_sampler << " " << cluster->kd_blobs().size() << " " << shad_cluster.kd_blobs().size() << " " << niblobs << std::endl;
                }

                //     // FIXME: These two methods need to be added to the Cluster Facade.
                //     // They should set/get "cluster_id" from the "cluster_scalar" PC.
                //     // FIXME: above we add cluster_id to the "cluster_scalar" PC.  Do we
                //     // want to also try to copy over the "light" index entry?

                auto cc2 = original.merge(splits,orig_cluster);
                for (const auto& val : cc2) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
                orig_cluster->put_pcarray(cc2, "isolated", "perblob");

                auto cc3 = shadow.merge(shadow_splits,shadow_orig_cluster);
                for (const auto& val : cc3) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
                shadow_orig_cluster->put_pcarray(cc3, "isolated", "perblob");

            }

        // FIXME: do we need/want to copy over any PCs in the grouping?
        // Specifically optical light/flash/flashlight PCs?

        }
    }
}
