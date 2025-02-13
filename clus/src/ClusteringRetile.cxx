// Developers: see important comments in header file.
//
// A "FIXME" really must be fixed if you expect anything to work.  Lower case
// "fixme" are "normal" fixmes that everyone ignores.


#include "WireCellClus/ClusteringRetile.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointTree.h"

using namespace WireCell;

// Segregate this weird choice for namespace.
namespace WCC = WireCell::PointCloud::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;


WCC::ClusteringRetile::ClusteringRetile(const WireCell::Configuration& cfg)
{
    auto sampler = get<std::string>(cfg, "sampler","");
    if (sampler.empty()) {
        raise<ValueError>("ClusteringRetile requires an IBlobSampler type/name");
    }
    auto face = get<std::string>(cfg, "face","");
    if (face.empty()) {
        raise<ValueError>("ClusteringRetile requires an IAnodeFace type/name");
    }

    m_sampler = Factory::find_tn<IBlobSampler>(sampler);
    m_face = Factory::find_tn<IAnodeFace>(face);
}


// Step 1. Build activities from blobs in a cluster.
WRG::activities_t WCC::ClusteringRetile::get_activity(const Cluster& cluster) const
{
    const double hit=1.0;       // actual charge value does not matter to tiling.

    WRG::activities_t ret;

    // first the horizontal/vertical bounds
    ret.emplace_back(0, 1, hit);
    ret.emplace_back(1, 1, hit);

    // fixme: this assumes "iend" is the usual one-past-last aka forming thigh
    // side of a half-open range.  I'm not sure if PointTreeBuilding is
    // violating that or not.
    auto tedious = [&](int ilayer, int ibeg, int iend) { ret.emplace_back(ilayer, iend-ibeg, hit, ibeg); };

    for (const auto& fblob : cluster.children()) {
        tedious(2, fblob->u_wire_index_min(), fblob->u_wire_index_max());
        tedious(3, fblob->v_wire_index_min(), fblob->v_wire_index_max());
        tedious(4, fblob->w_wire_index_min(), fblob->w_wire_index_max());
    }

    return ret;
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

    // Do the actual tiling.
    auto bshapes = WRG::make_blobs(m_face->raygrid(), activity);

    // Convert RayGrid blob shapes into IBlobs 
    const float blob_value = 0.0;  // tiling doesn't consider particular charge
    const float blob_error = 0.0;  // tiling doesn't consider particular charge
    int blob_ident=0;
    for (const auto& bshape : bshapes) {
        ISlice::pointer slice = nullptr; // fixme: okay?
        IBlob::pointer iblob = std::make_shared<Aux::SimpleBlob>(blob_ident++, blob_value,
                                                                 blob_error, bshape, slice, m_face);
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

        auto& shad_cluster = shadow.make_child();

        // FIXME: These two methods need to be added to the Cluster Facade.
        // They should set/get "cluster_id" from the "cluster_scalar" PC.
        shad_cluster.set_ident(orig_cluster->ident());

        // Step 1.
        auto orig_activity = get_activity(*orig_cluster);

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
        }

        // FIXME: above we add cluster_id to the "cluster_scalar" PC.  Do we
        // want to also try to copy over the "light" index entry?
    }
    // FIXME: do we need/want to copy over any PCs in the grouping?
    // Specifically optical light/flash/flashlight PCs?

}
