#include "WireCellHio/HIO.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"

#include <hdf5.h>
#include <algorithm>
#include <numeric>

using namespace WireCell;

namespace WireCell::Hio {

    // Convert DataType enum to HDF5 native type
    hid_t dtype_to_hid(Hio::DataType dtype) {
        switch(dtype) {
            case Hio::DataType::int16:    return H5T_NATIVE_INT16;
            case Hio::DataType::int32:    return H5T_NATIVE_INT32;
            case Hio::DataType::int64:    return H5T_NATIVE_INT64;
            case Hio::DataType::float32:  return H5T_NATIVE_FLOAT;
            case Hio::DataType::float64:  return H5T_NATIVE_DOUBLE;
        }
        raise<IOError>("Unknown DataType enum value");
        return -1;  // unreachable
    }

    // Convert C++ type to HDF5 native type
    template<typename T> hid_t native_type();
    template<> hid_t native_type<int16_t>() { return H5T_NATIVE_INT16; }
    template<> hid_t native_type<int32_t>() { return H5T_NATIVE_INT32; }
    template<> hid_t native_type<int64_t>() { return H5T_NATIVE_INT64; }
    template<> hid_t native_type<float>()    { return H5T_NATIVE_FLOAT; }
    template<> hid_t native_type<double>()   { return H5T_NATIVE_DOUBLE; }

    // Check HDF5 return value and throw if error
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

    // Create intermediate groups in path if needed
    void ensure_parents(hid_t file_id, const std::string& datapath) {
        // Split path into components
        std::vector<std::string> parts;
        size_t pos = 0;

        // Skip leading slash
        if (!datapath.empty() && datapath[0] == '/') {
            pos = 1;
        }

        while (pos < datapath.size()) {
            size_t next = datapath.find('/', pos);
            if (next == std::string::npos) {
                // Last component is the dataset name, not a group
                break;
            }
            if (next > pos) {
                parts.push_back(datapath.substr(pos, next - pos));
            }
            pos = next + 1;
        }

        // Create each group in the path
        std::string current_path;
        for (const auto& part : parts) {
            current_path += "/" + part;

            // Check if group exists
            H5E_BEGIN_TRY {
                hid_t existing = H5Gopen(file_id, current_path.c_str(), H5P_DEFAULT);
                if (existing >= 0) {
                    H5Gclose(existing);
                    continue;
                }
            } H5E_END_TRY;

            // Create group
            hid_t group = H5Gcreate(file_id, current_path.c_str(),
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            check_h5(group, String::format("Failed to create group: %s", current_path));
            H5Gclose(group);
        }
    }

    void write_json_attribute(hid_t loc_id, const std::string& attr_name,
                             const Json::Value& value) {
        // Delete attribute if it exists
        H5E_BEGIN_TRY {
            H5Adelete(loc_id, attr_name.c_str());
        } H5E_END_TRY;

        if (value.isNull()) {
            return;  // Don't write null values
        }

        if (value.isBool()) {
            hid_t atype = H5Tcopy(H5T_NATIVE_HBOOL);
            hid_t space = H5Screate(H5S_SCALAR);
            hid_t attr = H5Acreate(loc_id, attr_name.c_str(), atype, space,
                                  H5P_DEFAULT, H5P_DEFAULT);
            check_h5(attr, String::format("Failed to create bool attribute: %s", attr_name));
            hbool_t val = value.asBool();
            check_herr(H5Awrite(attr, atype, &val),
                      String::format("Failed to write bool attribute: %s", attr_name));
            H5Aclose(attr);
            H5Sclose(space);
            H5Tclose(atype);
        }
        else if (value.isInt()) {
            hid_t space = H5Screate(H5S_SCALAR);
            hid_t attr = H5Acreate(loc_id, attr_name.c_str(), H5T_NATIVE_INT64,
                                  space, H5P_DEFAULT, H5P_DEFAULT);
            check_h5(attr, String::format("Failed to create int attribute: %s", attr_name));
            int64_t val = value.asInt64();
            check_herr(H5Awrite(attr, H5T_NATIVE_INT64, &val),
                      String::format("Failed to write int attribute: %s", attr_name));
            H5Aclose(attr);
            H5Sclose(space);
        }
        else if (value.isDouble()) {
            hid_t space = H5Screate(H5S_SCALAR);
            hid_t attr = H5Acreate(loc_id, attr_name.c_str(), H5T_NATIVE_DOUBLE,
                                  space, H5P_DEFAULT, H5P_DEFAULT);
            check_h5(attr, String::format("Failed to create double attribute: %s", attr_name));
            double val = value.asDouble();
            check_herr(H5Awrite(attr, H5T_NATIVE_DOUBLE, &val),
                      String::format("Failed to write double attribute: %s", attr_name));
            H5Aclose(attr);
            H5Sclose(space);
        }
        else if (value.isString()) {
            std::string str = value.asString();
            hid_t atype = H5Tcopy(H5T_C_S1);
            H5Tset_size(atype, str.size());
            H5Tset_strpad(atype, H5T_STR_NULLTERM);
            hid_t space = H5Screate(H5S_SCALAR);
            hid_t attr = H5Acreate(loc_id, attr_name.c_str(), atype, space,
                                  H5P_DEFAULT, H5P_DEFAULT);
            check_h5(attr, String::format("Failed to create string attribute: %s", attr_name));
            check_herr(H5Awrite(attr, atype, str.c_str()),
                      String::format("Failed to write string attribute: %s", attr_name));
            H5Aclose(attr);
            H5Sclose(space);
            H5Tclose(atype);
        }
        else if (value.isArray() || value.isObject()) {
            // For complex types, serialize to JSON string
            Json::StreamWriterBuilder builder;
            builder["indentation"] = "";
            std::string json_str = Json::writeString(builder, value);

            hid_t atype = H5Tcopy(H5T_C_S1);
            H5Tset_size(atype, json_str.size());
            H5Tset_strpad(atype, H5T_STR_NULLTERM);
            hid_t space = H5Screate(H5S_SCALAR);
            hid_t attr = H5Acreate(loc_id, attr_name.c_str(), atype, space,
                                  H5P_DEFAULT, H5P_DEFAULT);
            check_h5(attr, String::format("Failed to create JSON attribute: %s", attr_name));
            check_herr(H5Awrite(attr, atype, json_str.c_str()),
                      String::format("Failed to write JSON attribute: %s", attr_name));
            H5Aclose(attr);
            H5Sclose(space);
            H5Tclose(atype);
        }
    }

    // Helper to read JSON value from HDF5 attribute
    Json::Value read_json_attribute(hid_t loc_id, const std::string& attr_name) {
        hid_t attr = H5Aopen(loc_id, attr_name.c_str(), H5P_DEFAULT);
        if (attr < 0) {
            return Json::Value();  // Attribute doesn't exist
        }

        hid_t atype = H5Aget_type(attr);
        H5T_class_t type_class = H5Tget_class(atype);

        Json::Value result;

        if (type_class == H5T_INTEGER) {
            if (H5Tequal(atype, H5T_NATIVE_HBOOL)) {
                hbool_t val;
                H5Aread(attr, atype, &val);
                result = static_cast<bool>(val);
            } else {
                int64_t val;
                H5Aread(attr, H5T_NATIVE_INT64, &val);
                result = static_cast<Json::Int64>(val);
            }
        }
        else if (type_class == H5T_FLOAT) {
            double val;
            H5Aread(attr, H5T_NATIVE_DOUBLE, &val);
            result = val;
        }
        else if (type_class == H5T_STRING) {
            size_t size = H5Tget_size(atype);
            std::vector<char> buffer(size + 1);
            H5Aread(attr, atype, buffer.data());
            buffer[size] = '\0';
            std::string str(buffer.data());

            // Try to parse as JSON for complex types
            if (!str.empty() && (str[0] == '{' || str[0] == '[')) {
                Json::CharReaderBuilder builder;
                Json::CharReader* reader = builder.newCharReader();
                std::string errors;
                if (reader->parse(str.c_str(), str.c_str() + str.size(),
                                 &result, &errors)) {
                    delete reader;
                } else {
                    delete reader;
                    result = str;  // Not JSON, just a string
                }
            } else {
                result = str;
            }
        }

        H5Tclose(atype);
        H5Aclose(attr);
        return result;
    }

    void show_errors(bool on)
    {
        if (on) {
            H5Eset_auto(H5E_DEFAULT, (H5E_auto_t)H5Eprint2, stderr);
        }
        else {
            H5Eset_auto(H5E_DEFAULT, NULL, NULL);
        }
    }

    std::vector<hsize_t> compute_chunks(const std::vector<int>& chunks,
                                        const std::vector<int64_t>& shape,
                                        size_t element_size,
                                        size_t target_bytes)
    {
        if (shape.empty()) {
            raise<IOError>("Cannot compute chunks for empty shape");
        }

        // Work with a mutable copy, reversed
        std::vector<hsize_t> result;
        for (auto it = chunks.rbegin(); it != chunks.rend(); ++it) {
            result.push_back(static_cast<hsize_t>(*it));
        }

        // If chunks empty, push back the size of the last dimension
        if (result.empty()) {
            result.push_back(static_cast<hsize_t>(shape[shape.size() - 1]));
        }

        // If chunks has size 1, push back a number so that this number
        // multiplied by the value in chunks and element_size is less than target
        if (result.size() == 1) {
            size_t last_chunk_bytes = result[0] * element_size;
            if (last_chunk_bytes > 0 && shape.size() >= 2) {
                // Calculate how many of the second-to-last dimension we can fit
                size_t second_last_chunk = target_bytes / last_chunk_bytes;
                // Clamp to actual dimension size
                hsize_t dim_size = static_cast<hsize_t>(shape[shape.size() - 2]);
                if (second_last_chunk > dim_size) {
                    second_last_chunk = dim_size;
                }
                if (second_last_chunk < 1) {
                    second_last_chunk = 1;
                }
                result.push_back(static_cast<hsize_t>(second_last_chunk));
            }
        }

        // Push back 1 until result has length equal to number of dimensions
        while (result.size() < shape.size()) {
            result.push_back(1);
        }

        // Reverse again to get final chunk sizes
        std::reverse(result.begin(), result.end());

        // Ensure chunks don't exceed shape dimensions
        for (size_t i = 0; i < result.size(); ++i) {
            if (result[i] > static_cast<hsize_t>(shape[i])) {
                result[i] = static_cast<hsize_t>(shape[i]);
            }
            if (result[i] < 1) {
                result[i] = 1;
            }
        }

        return result;
    }


    hid_t open(const std::string& filename, unsigned mode) {
        hid_t file_id = -1;

        if (mode == FileMode::trunc || mode == FileMode::excl) {
            file_id = H5Fcreate(filename.c_str(), mode, H5P_DEFAULT, H5P_DEFAULT);
            check_h5(file_id, String::format("Failed to create HDF5 file: %s", filename));
        }
        else {
            file_id = H5Fopen(filename.c_str(), mode, H5P_DEFAULT);
            check_h5(file_id, String::format("Failed to open HDF5 file: %s", filename));
        }

        return file_id;
    }

    void close(hid_t file_id) {
        check_herr(H5Fclose(file_id), "Failed to close HDF5 file");
    }

    void write_dataset(hid_t file_id, const std::byte* data,
                      const std::vector<int64_t>& shape,
                      DataType dtype,
                      const std::string& datapath,
                      int gzip,
                      const std::vector<int>& chunks) {
        if (shape.empty()) {
            raise<IOError>("Dataset shape cannot be empty");
        }

        ensure_parents(file_id, datapath);

        // Convert shape to hsize_t
        std::vector<hsize_t> dims(shape.begin(), shape.end());

        // Create dataspace
        hid_t space = H5Screate_simple(dims.size(), dims.data(), nullptr);
        check_h5(space, "Failed to create dataspace");

        // Get HDF5 type
        hid_t h5type = dtype_to_hid(dtype);

        // Determine element size
        size_t element_size = 0;
        switch(dtype) {
            case DataType::int16:   element_size = 2; break;
            case DataType::int32:   element_size = 4; break;
            case DataType::int64:   element_size = 8; break;
            case DataType::float32: element_size = 4; break;
            case DataType::float64: element_size = 8; break;
        }

        // Setup dataset creation properties
        hid_t dcpl = H5P_DEFAULT;

        // Clamp gzip to valid range [0, 9]
        int gzip_level = gzip;
        if (gzip_level < 0) gzip_level = 0;
        if (gzip_level > 9) gzip_level = 9;

        if (gzip_level > 0) {
            // Compute chunk sizes
            std::vector<hsize_t> chunk_dims = compute_chunks(chunks, shape, element_size);

            // Create dataset creation property list
            dcpl = H5Pcreate(H5P_DATASET_CREATE);
            check_h5(dcpl, "Failed to create dataset creation property list");

            // Set chunk sizes
            check_herr(H5Pset_chunk(dcpl, chunk_dims.size(), chunk_dims.data()),
                      "Failed to set chunk sizes");

            // Set deflate (gzip) compression
            check_herr(H5Pset_deflate(dcpl, gzip_level),
                      "Failed to set deflate compression");
        }

        // Create dataset
        hid_t dset = H5Dcreate(file_id, datapath.c_str(), h5type, space,
                              H5P_DEFAULT, dcpl, H5P_DEFAULT);
        check_h5(dset, String::format("Failed to create dataset: %s", datapath));

        // Write data
        check_herr(H5Dwrite(dset, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
                  String::format("Failed to write dataset: %s", datapath));

        H5Dclose(dset);
        H5Sclose(space);

        if (dcpl != H5P_DEFAULT) {
            H5Pclose(dcpl);
        }
    }

    void read_dataset(hid_t file_id, std::vector<std::byte>& data,
                     std::vector<int64_t>& shape,
                     const std::string& datapath) {
        // Open dataset
        hid_t dset = H5Dopen(file_id, datapath.c_str(), H5P_DEFAULT);
        check_h5(dset, String::format("Failed to open dataset: %s", datapath));

        // Get dataspace
        hid_t space = H5Dget_space(dset);
        check_h5(space, String::format("Failed to get dataspace for: %s", datapath));

        // Get dimensions
        int ndims = H5Sget_simple_extent_ndims(space);
        std::vector<hsize_t> dims(ndims);
        H5Sget_simple_extent_dims(space, dims.data(), nullptr);
        shape.assign(dims.begin(), dims.end());

        // Get type and size
        hid_t dtype = H5Dget_type(dset);
        size_t type_size = H5Tget_size(dtype);

        // Calculate total size
        size_t total_elements = std::accumulate(dims.begin(), dims.end(),
                                               size_t(1), std::multiplies<size_t>());
        size_t total_bytes = total_elements * type_size;

        // Resize data buffer
        data.resize(total_bytes);

        // Read data
        check_herr(H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()),
                  String::format("Failed to read dataset: %s", datapath));

        H5Tclose(dtype);
        H5Sclose(space);
        H5Dclose(dset);
    }

    void write_metadata(hid_t file_id, const Configuration& metadata,
                       const std::string& datapath) {
        if (metadata.isNull() || metadata.empty()) {
            return;  // Nothing to write
        }

        // Try to open as dataset first, then as group
        hid_t loc_id = -1;
        bool is_dataset = false;

        H5E_BEGIN_TRY {
            loc_id = H5Dopen(file_id, datapath.c_str(), H5P_DEFAULT);
            if (loc_id >= 0) {
                is_dataset = true;
            }
        } H5E_END_TRY;

        if (loc_id < 0) {
            H5E_BEGIN_TRY {
                loc_id = H5Gopen(file_id, datapath.c_str(), H5P_DEFAULT);
            } H5E_END_TRY;
        }

        check_h5(loc_id, String::format("Failed to open object for metadata: %s", datapath));

        // Write each metadata field as an attribute
        for (const auto& key : metadata.getMemberNames()) {
            write_json_attribute(loc_id, key, metadata[key]);
        }

        if (is_dataset) {
            H5Dclose(loc_id);
        } else {
            H5Gclose(loc_id);
        }
    }

    Configuration read_metadata(hid_t file_id, const std::string& datapath) {
        // Try to open as dataset first, then as group
        hid_t loc_id = -1;
        bool is_dataset = false;

        H5E_BEGIN_TRY {
            loc_id = H5Dopen(file_id, datapath.c_str(), H5P_DEFAULT);
            if (loc_id >= 0) {
                is_dataset = true;
            }
        } H5E_END_TRY;

        if (loc_id < 0) {
            H5E_BEGIN_TRY {
                loc_id = H5Gopen(file_id, datapath.c_str(), H5P_DEFAULT);
            } H5E_END_TRY;
        }

        if (loc_id < 0) {
            return Configuration();  // Object doesn't exist, return empty
        }

        Configuration result = Json::objectValue;

        // Get number of attributes
#if H5_VERSION_GE(1, 12, 0)
        H5O_info2_t obj_info;
        check_herr(H5Oget_info3(loc_id, &obj_info, H5O_INFO_NUM_ATTRS),
                  String::format("Failed to get object info: %s", datapath));
#else
        H5O_info_t obj_info;
        check_herr(H5Oget_info(loc_id, &obj_info),
                  String::format("Failed to get object info: %s", datapath));
#endif

        // Read each attribute
        for (hsize_t i = 0; i < obj_info.num_attrs; ++i) {
            char attr_name[256];
            H5Aget_name_by_idx(loc_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i,
                              attr_name, sizeof(attr_name), H5P_DEFAULT);
            result[attr_name] = read_json_attribute(loc_id, attr_name);
        }

        if (is_dataset) {
            H5Dclose(loc_id);
        } else {
            H5Gclose(loc_id);
        }

        return result;
    }

    // Template explicit instantiations for write_dataset
    template<typename T>
    void write_dataset(hid_t file_id, const std::vector<T>& data,
                      const std::vector<int64_t>& shape,
                      const std::string& datapath,
                      int gzip,
                      const std::vector<int>& chunks) {
        if (shape.empty()) {
            raise<IOError>("Dataset shape cannot be empty");
        }

        // Verify shape matches data size
        size_t expected_size = std::accumulate(shape.begin(), shape.end(),
                                              size_t(1), std::multiplies<size_t>());
        if (data.size() != expected_size) {
            raise<IOError>(String::format(
                "Data size %zu does not match shape size %zu for dataset: %s",
                data.size(), expected_size, datapath));
        }

        ensure_parents(file_id, datapath);

        // Convert shape to hsize_t
        std::vector<hsize_t> dims(shape.begin(), shape.end());

        // Create dataspace
        hid_t space = H5Screate_simple(dims.size(), dims.data(), nullptr);
        check_h5(space, "Failed to create dataspace");

        // Get HDF5 type
        hid_t h5type = native_type<T>();

        // Setup dataset creation properties
        hid_t dcpl = H5P_DEFAULT;

        // Clamp gzip to valid range [0, 9]
        int gzip_level = gzip;
        if (gzip_level < 0) gzip_level = 0;
        if (gzip_level > 9) gzip_level = 9;

        if (gzip_level > 0) {
            // Compute chunk sizes
            size_t element_size = sizeof(T);
            std::vector<hsize_t> chunk_dims = compute_chunks(chunks, shape, element_size);

            // Create dataset creation property list
            dcpl = H5Pcreate(H5P_DATASET_CREATE);
            check_h5(dcpl, "Failed to create dataset creation property list");

            // Set chunk sizes
            check_herr(H5Pset_chunk(dcpl, chunk_dims.size(), chunk_dims.data()),
                      "Failed to set chunk sizes");

            // Set deflate (gzip) compression
            check_herr(H5Pset_deflate(dcpl, gzip_level),
                      "Failed to set deflate compression");
        }

        // Create dataset
        hid_t dset = H5Dcreate(file_id, datapath.c_str(), h5type, space,
                              H5P_DEFAULT, dcpl, H5P_DEFAULT);
        check_h5(dset, String::format("Failed to create dataset: %s", datapath));

        // Write data
        check_herr(H5Dwrite(dset, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()),
                  String::format("Failed to write dataset: %s", datapath));

        H5Dclose(dset);
        H5Sclose(space);

        if (dcpl != H5P_DEFAULT) {
            H5Pclose(dcpl);
        }
    }

    // Template explicit instantiations for read_dataset
    template<typename T>
    void read_dataset(hid_t file_id, std::vector<T>& data,
                     std::vector<int64_t>& shape,
                     const std::string& datapath) {
        // Open dataset
        hid_t dset = H5Dopen(file_id, datapath.c_str(), H5P_DEFAULT);
        check_h5(dset, String::format("Failed to open dataset: %s", datapath));

        // Get dataspace
        hid_t space = H5Dget_space(dset);
        check_h5(space, String::format("Failed to get dataspace for: %s", datapath));

        // Get dimensions
        int ndims = H5Sget_simple_extent_ndims(space);
        std::vector<hsize_t> dims(ndims);
        H5Sget_simple_extent_dims(space, dims.data(), nullptr);
        shape.assign(dims.begin(), dims.end());

        // Calculate total size
        size_t total_elements = std::accumulate(dims.begin(), dims.end(),
                                               size_t(1), std::multiplies<size_t>());

        // Resize data buffer
        data.resize(total_elements);

        // Get HDF5 type
        hid_t h5type = native_type<T>();

        // Read data
        check_herr(H5Dread(dset, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()),
                  String::format("Failed to read dataset: %s", datapath));

        H5Sclose(space);
        H5Dclose(dset);
    }

    // Explicit instantiations
    template void write_dataset<int16_t>(hid_t, const std::vector<int16_t>&,
                                        const std::vector<int64_t>&,
                                        const std::string&,
                                        int, const std::vector<int>&);
    template void write_dataset<int32_t>(hid_t, const std::vector<int32_t>&,
                                        const std::vector<int64_t>&,
                                        const std::string&,
                                        int, const std::vector<int>&);
    template void write_dataset<int64_t>(hid_t, const std::vector<int64_t>&,
                                        const std::vector<int64_t>&,
                                        const std::string&,
                                        int, const std::vector<int>&);
    template void write_dataset<float>(hid_t, const std::vector<float>&,
                                      const std::vector<int64_t>&,
                                      const std::string&,
                                      int, const std::vector<int>&);
    template void write_dataset<double>(hid_t, const std::vector<double>&,
                                       const std::vector<int64_t>&,
                                       const std::string&,
                                       int, const std::vector<int>&);

    template void read_dataset<int16_t>(hid_t, std::vector<int16_t>&,
                                       std::vector<int64_t>&,
                                       const std::string&);
    template void read_dataset<int32_t>(hid_t, std::vector<int32_t>&,
                                       std::vector<int64_t>&,
                                       const std::string&);
    template void read_dataset<int64_t>(hid_t, std::vector<int64_t>&,
                                       std::vector<int64_t>&,
                                       const std::string&);
    template void read_dataset<float>(hid_t, std::vector<float>&,
                                     std::vector<int64_t>&,
                                     const std::string&);
    template void read_dataset<double>(hid_t, std::vector<double>&,
                                      std::vector<int64_t>&,
                                      const std::string&);

}  // namespace WireCell::Hio
