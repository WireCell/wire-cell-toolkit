/**
   Utilities support the tensor data model for point cloud trees.
 */
#ifndef WIRECELLAUX_TENSORDM
#define WIRECELLAUX_TENSORDM

#include "WireCellAux/TensorDMdataset.h" // as_tensors
#include "WireCellAux/TensorDMcommon.h" // make_metadata_tensor
#include "WireCellUtil/PointTree.h"
#include "WireCellIface/ITensor.h"

namespace WireCell::Aux::TensorDM {

    /// Named point cloud sets 

    /** Convert a named point cloud set to tensors.

        The first tensor in the returned vector will represent the
        pcnamedset and subsequent entries will represent each item in
        the set.

        If <store> is not given it will be set to <datapath>/namedpcs.
        The pcdataset tensors will have datapaths of <store>/<name>.

     */
    template<typename PairIter>
    ITensor::vector pcnamedset_as_tensors(PairIter beg, PairIter end,
                                          const std::string& datapath,
                                          std::string store = "") {
        if (store.empty()) {
            store = datapath + "/namedpcs";
        }

        Configuration md, items = Json::objectValue;

        ITensor::vector ret;
        ret.push_back(nullptr);     // fill in below
        for (auto it = beg; it != end; ++it) {
            const std::string& name = it->first;
            const auto& ds = it->second;
            auto dp = store + "/" + name;
            items[name] = dp;
            auto tens = as_tensors(ds, dp);
            ret.insert(ret.end(), tens.begin(), tens.end());
        }
        md["items"] = items;
        ret[0] = make_metadata_tensor("pcnamedset", datapath, md);
        return ret;
    }

    using named_pointclouds_t = WireCell::PointCloud::Tree::named_pointclouds_t;

    inline
    ITensor::vector as_tensors(const named_pointclouds_t& pcns,
                               const std::string& datapath,
                               std::string store = "") {
        return pcnamedset_as_tensors(pcns.begin(), pcns.end(), datapath, store);
    }

    
    /** Convert tensors representing a named point cloud set.

       The tensor found at datapath is converted unless datapath is
       empty in which case the first pcnamedset found is used.

     */
    named_pointclouds_t as_pcnamedset(const ITensor::vector& tens, const std::string& datapath = "");
    named_pointclouds_t as_pcnamedset(const TensorIndex& ti, const std::string& datapath = "");
                                                             


    /**
     * Convert a point cloud tree to vector of ITensors.
     */ 
    ITensor::vector as_tensors(
        const WireCell::PointCloud::Tree::Points::node_t& root,
        const std::string& datapath);

    /**
     * Convert sequence of ITensor to a point cloud tree.
     *
     * The datapath specifies which in the sequence of tensors is the
     * pctree.  If empty then the first tensor of type pctree is
     * assumed.
     */
    std::unique_ptr<WireCell::PointCloud::Tree::Points::node_t>
    as_pctree(const ITensor::vector& tens,
              const std::string& datapath="");
    std::unique_ptr<WireCell::PointCloud::Tree::Points::node_t>
    as_pctree(const TensorIndex& ti,
              const std::string& datapath="");

}

#endif
