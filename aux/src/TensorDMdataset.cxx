#include "WireCellAux/TensorDM.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Type.h"
#include <boost/filesystem.hpp>

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;


template<typename ElementType>
ITensor::pointer make_tensor(const PointCloud::Array& array, const std::string& datapath)
{
    ElementType* data = (ElementType*)array.bytes().data();
    auto md = array.metadata();
    md["datapath"] = datapath;
    md["datatype"] = "pcarray";
    return std::make_shared<SimpleTensor>(array.shape(), data, md);
}

ITensor::pointer WireCell::Aux::TensorDM::as_tensor(const PointCloud::Array& array,
                                                    const std::string& datapath)
{
    if (array.is_type<char>()) return make_tensor<char>(array, datapath);
    if (array.is_type<std::byte>()) return make_tensor<std::byte>(array, datapath);

    if (array.is_type<int8_t>())  return make_tensor<int8_t>(array, datapath);
    if (array.is_type<int16_t>()) return make_tensor<int16_t>(array, datapath);
    if (array.is_type<int32_t>()) return make_tensor<int32_t>(array, datapath);
    if (array.is_type<int64_t>()) return make_tensor<int64_t>(array, datapath);

    if (array.is_type<uint8_t>()) return make_tensor<uint8_t>(array, datapath);
    if (array.is_type<uint16_t>()) return make_tensor<uint16_t>(array, datapath);
    if (array.is_type<uint32_t>()) return make_tensor<uint32_t>(array, datapath);
    if (array.is_type<uint64_t>()) return make_tensor<uint64_t>(array, datapath);

    if (array.is_type<float>()) return make_tensor<float>(array, datapath);
    if (array.is_type<double>()) return make_tensor<double>(array, datapath);
    if (array.is_type<std::complex<float>>()) return make_tensor<std::complex<float>>(array, datapath);
    if (array.is_type<std::complex<double>>()) return make_tensor<std::complex<double>>(array, datapath);

    THROW(ValueError() << errmsg{"unsupported point cloud array type: " + array.dtype()});
}

ITensor::vector
WireCell::Aux::TensorDM::as_tensors(const PointCloud::Dataset& dataset,
                                    const std::string& datapath)
{
    auto md = dataset.metadata();
    md["datapath"] = datapath;
    md["datatype"] = "pcdataset";
    md["arrays"] = Json::objectValue;

    boost::filesystem::path trpath = datapath + "/arrays";
    ITensor::vector ret;

    int arrnum = 0;
    for (const auto& [name, arr] : dataset.store()) {
        auto dpath = trpath / std::to_string(arrnum++);
        ret.push_back(as_tensor(arr, dpath.string()));
        md["arrays"][name] = dpath.string();
    }

    auto top = std::make_shared<SimpleTensor>(md);
    ret.insert(ret.begin(), top);
    return ret;
}

template<typename ElementType>
PointCloud::Array make_array(const ITensor::pointer& ten, bool share)
{
    ElementType* data = (ElementType*)ten->data();
    PointCloud::Array arr(data, ten->shape(), share);
    arr.metadata() = ten->metadata();
    return arr;
}

PointCloud::Array WireCell::Aux::TensorDM::as_array(const ITensor::pointer& ten, bool share)
{
    const auto& dt = ten->element_type();

    if (dt == typeid(char)) return make_array<char>(ten, share);
    if (dt == typeid(std::byte)) return make_array<std::byte>(ten, share);

    if (dt == typeid(int8_t)) return make_array<int8_t>(ten, share);
    if (dt == typeid(int16_t)) return make_array<int16_t>(ten, share);
    if (dt == typeid(int32_t)) return make_array<int32_t>(ten, share);
    if (dt == typeid(int64_t)) return make_array<int64_t>(ten, share);

    if (dt == typeid(uint8_t)) return make_array<uint8_t>(ten, share);
    if (dt == typeid(uint16_t)) return make_array<uint16_t>(ten, share);
    if (dt == typeid(uint32_t)) return make_array<uint32_t>(ten, share);
    if (dt == typeid(uint64_t)) return make_array<uint64_t>(ten, share);

    if (dt == typeid(float)) return make_array<float>(ten, share);            
    if (dt == typeid(double)) return make_array<double>(ten, share);            
    if (dt == typeid(std::complex<float>)) return make_array<std::complex<float>>(ten, share);            
    if (dt == typeid(std::complex<double>)) return make_array<std::complex<double>>(ten, share);            

    THROW(ValueError() << errmsg{"unsupported ITensor type " + demangle(dt.name())});
}

PointCloud::Dataset
WireCell::Aux::TensorDM::as_dataset(const ITensor::vector& tens,
                                    const std::string datapath,
                                    bool share)
{

    ITensor::pointer top = nullptr;
    std::unordered_map<std::string, ITensor::pointer> located;
    for (const auto& iten : tens) {
        const auto& tenmd = iten->metadata();
        const auto dtype = tenmd["datatype"].asString();
        const auto dpath = tenmd["datapath"].asString();
        if (!top and dtype == "pcdataset") {
            if (datapath.empty() or datapath == dpath) {
                top = iten;
            }
            continue;
        }
        if (dtype == "pcarray") {
            located[dpath] = iten;
        }
        continue;
    }

    if (!top) {
        THROW(ValueError() << errmsg{"no array of datatype \"pcdataset\""});
    }

    PointCloud::Dataset ret;
    auto topmd = top->metadata();
    ret.metadata() = topmd;
    auto arrrefs = topmd["arrays"];
    for (const auto& name : arrrefs.getMemberNames()) {
        const auto path = arrrefs[name].asString();
        auto it = located.find(path);
        if (it == located.end()) {
            THROW(ValueError() << errmsg{"no array \"" + name + "\" at path \"" + path + "\""});
        }
        auto arr = as_array(it->second, share);
        ret.add(name, arr);
    }

    return ret;    
}


PointCloud::Dataset
WireCell::Aux::TensorDM::as_dataset(const ITensorSet::pointer& tens,
                                    const std::string datapath,
                                    bool share)
{
    auto sv = tens->tensors();
    return as_dataset(*(sv.get()), datapath, share);
}
