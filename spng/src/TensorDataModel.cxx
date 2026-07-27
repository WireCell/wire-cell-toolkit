#include "WireCellSpng/TensorDataModel.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include <regex>

namespace WireCell::SPNG {


    ITorchTensor::vector find_tensors(const ITorchTensorSet::pointer& ts, const std::string& datapath_regex)
    {
        auto match = std::regex(datapath_regex);
        ITorchTensor::vector found;
        for (auto iten : *(ts->tensors())) {
            const auto dp = datapath(iten);
            if (std::regex_match(dp, match)) {
                found.push_back(iten);
            }
        }
        return found;
    }


    ITorchTensor::pointer make_tensor(const std::string& datatype,
                                      const std::string& datapath,
                                      Configuration metadata,
                                      torch::Tensor ten,
                                      const std::string& parent,
                                      size_t batches)
    {
        metadata["datatype"] = datatype;
        metadata["datapath"] = datapath;
        if (! parent.empty()) {
            metadata["parent"] = parent;
        }
        if (batches) {
            metadata["batches"] = batches;
        }
        return std::make_shared<SimpleTorchTensor>(ten, metadata);
    }

    ITorchTensor::pointer make_frame(const std::string& base_datapath, int ident, double time, double period)
    {
        Configuration md;
        md["ident"] = ident;
        md["time"] = time;
        md["period"] = period;
        return make_tensor("frame", base_datapath + "/frame", md);
    }


    ITorchTensor::pointer make_traces(const std::string& base_datapath,
                                      const std::string& tag,
                                      torch::Tensor ten,
                                      int tbin)
    {
        Configuration md;
        md["tag"] = tag;
        md["tbin"] = tbin;
        int ndim = ten.dim();
        size_t batches = 0;
        if (ndim > 2) {
            batches = ten.size(0);
        }
        return make_tensor("traces", base_datapath + "/traces/" + tag,
                           md, ten, base_datapath + "/frame", batches);
    }


    ITorchTensor::pointer make_chids(const std::string& base_datapath,
                                     const std::string& tag,
                                     torch::Tensor ten)
    {
        Configuration md;
        md["tag"] = tag;
        int ndim = ten.dim();
        size_t batches = 0;
        if (ndim > 1) {
            batches = ten.size(0);
        }
        return make_tensor("chids", base_datapath + "/chids/" + tag,
                           md, ten, base_datapath + "/frame", batches);
        

    }

    /// Make a summaries tensor
    ITorchTensor::pointer make_summaries(const std::string& base_datapath,
                                         const std::string& tag,
                                         torch::Tensor ten)
    {
        Configuration md;
        md["tag"] = tag;
        int ndim = ten.dim();
        size_t batches = 0;
        if (ndim > 1) {
            batches = ten.size(0);
        }
        return make_tensor("summaries", base_datapath + "/summaries/" + tag,
                           md, ten, base_datapath + "/frame", batches);
    }
    
    /// Make a channel masks tensor
    ITorchTensor::pointer make_chmasks(const std::string& base_datapath,
                                       const std::string& label,
                                       torch::Tensor ten)
    {
        Configuration md;
        md["label"] = label;

        return make_tensor("chmasks", base_datapath + "/chmasks/" + label,
                           md, ten, base_datapath + "/frame");

    }

    ITorchTensor::vector make_frame_set(const std::string& base_datapath, IFrame::pointer iframe)
    {
        ITorchTensor::vector ret;
        return ret;             // not yet
    }


}
