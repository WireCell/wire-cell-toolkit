#include "WireCellHio/Tensors.h"
#include "WireCellHio/HIO.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"

#include <hdf5.h>

using namespace WireCell;

namespace WireCell::Hio {

    void check_h5(hid_t ret, const std::string& msg) {
        if (ret < 0) {
            raise<IOError>(msg);
        }
    }

    void check_herr(herr_t ret, const std::string& msg) {
        if (ret < 0) {
            raise<IOError>(msg);
        }
    }

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

}  // namespace WireCell::Hio
