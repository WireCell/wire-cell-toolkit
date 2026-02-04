#include "WireCellHio/Tensors.h"
#include "WireCellHio/HIO.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"

#include <hdf5.h>
#include <vector>

using namespace WireCell;

namespace WireCell::Hio {

    void write_itensor(hid_t file_id, ITensor::pointer tensor, const std::string& datapath) {
        if (!tensor) {
            raise<IOError>("Cannot write null ITensor pointer");
        }

        // Get tensor properties
        auto shape = tensor->shape();
        if (shape.empty()) {
            raise<IOError>(String::format("ITensor has empty shape for datapath: %s", datapath));
        }

        // Convert shape from size_t to int64_t
        std::vector<int64_t> h5_shape(shape.begin(), shape.end());

        // Get the raw data pointer
        const std::byte* data_ptr = tensor->data();
        size_t data_size = tensor->size();

        if (data_size == 0) {
            raise<IOError>(String::format("ITensor has zero data size for datapath: %s", datapath));
        }

        // Get dtype string and determine HDF5 type
        std::string dtype_str = tensor->dtype();
        if (dtype_str.empty()) {
            raise<IOError>(String::format("ITensor has unknown dtype for datapath: %s", datapath));
        }

        // Map dtype string to HDF5 DataType enum
        DataType dtype;
        if (dtype_str == "i2") {
            dtype = DataType::int16;
        } else if (dtype_str == "i4") {
            dtype = DataType::int32;
        } else if (dtype_str == "i8") {
            dtype = DataType::int64;
        } else if (dtype_str == "f4") {
            dtype = DataType::float32;
        } else if (dtype_str == "f8") {
            dtype = DataType::float64;
        } else {
            raise<IOError>(String::format(
                "Unsupported ITensor dtype '%s' for HDF5 serialization at datapath: %s",
                dtype_str, datapath));
        }

        // Write the dataset
        write_dataset(file_id, data_ptr, h5_shape, dtype, datapath);

        // Write metadata as attributes
        Configuration metadata = tensor->metadata();
        if (!metadata.isNull() && !metadata.empty()) {
            write_metadata(file_id, metadata, datapath);
        }
    }

    ITensor::pointer read_itensor(hid_t file_id, const std::string& datapath) {
        // Open dataset to query type
        hid_t dset = H5Dopen(file_id, datapath.c_str(), H5P_DEFAULT);
        check_h5(dset, String::format("Failed to open dataset: %s", datapath));

        // Get datatype
        hid_t dtype = H5Dget_type(dset);
        check_h5(dtype, String::format("Failed to get datatype for: %s", datapath));

        // Get dataspace and shape
        hid_t space = H5Dget_space(dset);
        check_h5(space, String::format("Failed to get dataspace for: %s", datapath));

        int ndims = H5Sget_simple_extent_ndims(space);
        std::vector<hsize_t> h5_dims(ndims);
        H5Sget_simple_extent_dims(space, h5_dims.data(), nullptr);

        ITensor::shape_t shape(h5_dims.begin(), h5_dims.end());

        // Read metadata
        Configuration metadata = read_metadata(file_id, datapath);

        // Determine the type and read accordingly
        ITensor::pointer result;
        H5T_class_t type_class = H5Tget_class(dtype);
        size_t type_size = H5Tget_size(dtype);

        if (type_class == H5T_INTEGER) {
            if (type_size == 2) {
                // int16
                std::vector<int16_t> data;
                std::vector<int64_t> tmp_shape;
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                read_dataset(file_id, data, tmp_shape, datapath);
                result = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);
            }
            else if (type_size == 4) {
                // int32
                std::vector<int32_t> data;
                std::vector<int64_t> tmp_shape;
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                read_dataset(file_id, data, tmp_shape, datapath);
                result = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);
            }
            else if (type_size == 8) {
                // int64
                std::vector<int64_t> data;
                std::vector<int64_t> tmp_shape;
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                read_dataset(file_id, data, tmp_shape, datapath);
                result = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);
            }
            else {
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                raise<IOError>(String::format(
                    "Unsupported integer size %zu bytes for datapath: %s",
                    type_size, datapath));
            }
        }
        else if (type_class == H5T_FLOAT) {
            if (type_size == 4) {
                // float32
                std::vector<float> data;
                std::vector<int64_t> tmp_shape;
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                read_dataset(file_id, data, tmp_shape, datapath);
                result = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);
            }
            else if (type_size == 8) {
                // float64
                std::vector<double> data;
                std::vector<int64_t> tmp_shape;
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                read_dataset(file_id, data, tmp_shape, datapath);
                result = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);
            }
            else {
                H5Sclose(space);
                H5Tclose(dtype);
                H5Dclose(dset);
                raise<IOError>(String::format(
                    "Unsupported float size %zu bytes for datapath: %s",
                    type_size, datapath));
            }
        }
        else {
            H5Sclose(space);
            H5Tclose(dtype);
            H5Dclose(dset);
            raise<IOError>(String::format(
                "Unsupported HDF5 type class for datapath: %s", datapath));
        }

        return result;
    }

    void write_itensorset(hid_t file_id, ITensorSet::pointer tensorset, const std::string& datapath) {
        if (!tensorset) {
            raise<IOError>("Cannot write null ITensorSet pointer");
        }

        ensure_parents(file_id, datapath);
        // // Ensure parent groups exist by extracting parent path and creating groups
        // std::string parent_path;
        // size_t last_slash = datapath.rfind('/');
        // if (last_slash != std::string::npos && last_slash > 0) {
        //     parent_path = datapath.substr(0, last_slash);

        //     // Create all parent groups
        //     std::vector<std::string> parts;
        //     size_t pos = 1; // Skip leading slash
        //     while (pos < parent_path.size()) {
        //         size_t next = parent_path.find('/', pos);
        //         if (next == std::string::npos) {
        //             parts.push_back(parent_path.substr(pos));
        //             break;
        //         }
        //         if (next > pos) {
        //             parts.push_back(parent_path.substr(pos, next - pos));
        //         }
        //         pos = next + 1;
        //     }

        //     // Create each parent group if it doesn't exist
        //     std::string current_path;
        //     for (const auto& part : parts) {
        //         current_path += "/" + part;

        //         // Check if group exists
        //         H5E_BEGIN_TRY {
        //             hid_t existing = H5Gopen(file_id, current_path.c_str(), H5P_DEFAULT);
        //             if (existing >= 0) {
        //                 H5Gclose(existing);
        //                 continue;
        //             }
        //         } H5E_END_TRY;

        //         // Create group
        //         hid_t parent_group = H5Gcreate(file_id, current_path.c_str(),
        //                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //         check_h5(parent_group, String::format("Failed to create parent group: %s", current_path));
        //         H5Gclose(parent_group);
        //     }
        // }

        // Create the group for the tensor set
        hid_t group_id = H5Gcreate(file_id, datapath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check_h5(group_id, String::format("Failed to create group: %s", datapath));

        // Write metadata as attributes on the group
        Configuration metadata = tensorset->metadata();
        if (!metadata.isNull() && !metadata.empty()) {
            write_metadata(file_id, metadata, datapath);
        }

        // Write ident as an attribute
        Configuration ident_config;
        ident_config["ident"] = tensorset->ident();
        write_metadata(file_id, ident_config, datapath);

        // Get tensors from the set
        auto tensors_ptr = tensorset->tensors();
        if (!tensors_ptr) {
            H5Gclose(group_id);
            raise<IOError>(String::format("ITensorSet has null tensors pointer at datapath: %s", datapath));
        }

        const auto& tensors = *tensors_ptr;

        // Write each tensor as a dataset in the group
        for (size_t i = 0; i < tensors.size(); ++i) {
            std::string tensor_path = datapath + "/" + std::to_string(i);
            write_itensor(file_id, tensors[i], tensor_path);
        }

        H5Gclose(group_id);
    }

    ITensorSet::pointer read_itensorset(hid_t file_id, const std::string& datapath) {
        // Open the group
        hid_t group_id = H5Gopen(file_id, datapath.c_str(), H5P_DEFAULT);
        check_h5(group_id, String::format("Failed to open group: %s", datapath));

        // Read metadata from the group
        Configuration metadata = read_metadata(file_id, datapath);

        // Extract ident from metadata
        int ident = 0;
        if (metadata.isMember("ident")) {
            ident = metadata["ident"].asInt();
            // Remove ident from metadata since it's a separate field
            metadata.removeMember("ident");
        }

        // Get number of objects in the group
#if H5_VERSION_GE(1, 12, 0)
        H5G_info_t group_info;
        check_herr(H5Gget_info(group_id, &group_info),
                  String::format("Failed to get group info: %s", datapath));
        hsize_t num_objs = group_info.nlinks;
#else
        hsize_t num_objs;
        check_herr(H5Gget_num_objs(group_id, &num_objs),
                  String::format("Failed to get number of objects: %s", datapath));
#endif

        // Read each tensor from the group
        auto tensors = std::make_shared<ITensor::vector>();
        for (hsize_t i = 0; i < num_objs; ++i) {
            // Get object name
#if H5_VERSION_GE(1, 12, 0)
            char name[256];
            ssize_t name_size = H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE,
                                                   i, name, sizeof(name), H5P_DEFAULT);
            check_h5(name_size, String::format("Failed to get object name at index %zu: %s", i, datapath));
#else
            char name[256];
            ssize_t name_size = H5Gget_objname_by_idx(group_id, i, name, sizeof(name));
            check_h5(name_size, String::format("Failed to get object name at index %zu: %s", i, datapath));
#endif

            // Read the tensor
            std::string tensor_path = datapath + "/" + std::string(name);
            ITensor::pointer tensor = read_itensor(file_id, tensor_path);
            tensors->push_back(tensor);
        }

        H5Gclose(group_id);

        // Construct and return SimpleTensorSet
        return std::make_shared<Aux::SimpleTensorSet>(ident, metadata, tensors);
    }

}  // namespace WireCell::Hio
