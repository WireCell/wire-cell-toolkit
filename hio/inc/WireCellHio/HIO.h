/** Thin wrapper around HDF5 C library
 *
 * HIO simplifies the hdf5.h API with minimal C++ wrapper.
 * Not all HDF5 functionality is exposed.
 *
 * Lives in the C++ namespace WireCell::Hio.
 *
 * The HIO API is based on functions, no objects other than STL containers.
 * The HDF5 hid_t opaque HDF5 ID is exposed by HIO functions.
 *
 * HIO functions will interpret any hdf5.h function return and will throw
 * the Wire-Cell Toolkit exception IOError with meaningful messages.
 */

#ifndef WIRECELLHIO_HIO
#define WIRECELLHIO_HIO

#include "WireCellUtil/Configuration.h"
#include <hdf5.h>

#include <cstddef>
#include <string>
#include <vector>
#include <cstdint>

namespace WireCell::Hio {

    /// File access modes
    struct FileMode {
        static inline const unsigned trunc  = H5F_ACC_TRUNC;   ///< Truncate existing file, create if doesn't exist
        static inline const unsigned excl   = H5F_ACC_EXCL;    ///< Fail if file exists
        static inline const unsigned rdrw   = H5F_ACC_RDWR;    ///< Read-write access to existing file
        static inline const unsigned rdonly = H5F_ACC_RDONLY;  ///< Read-only access to existing file
    };

    /// HDF5 data types
    enum struct DataType {
        int16,
        int32,
        int64,
        float32,
        float64
    };

    /// Convert DataType enum to HDF5 native type
    hid_t dtype_to_hid(Hio::DataType dtype);

    /** Check HDF5 return value and throw IOError if error
     *
     * @param ret HDF5 return value (hid_t)
     * @param msg Error message to include in exception
     * @throws IOError if ret < 0
     */
    void check_h5(hid_t ret, const std::string& msg);

    /** Check HDF5 error return and throw IOError if error
     *
     * @param ret HDF5 error return value (herr_t)
     * @param msg Error message to include in exception
     * @throws IOError if ret < 0
     */
    void check_herr(herr_t ret, const std::string& msg);

    /// Call to with false to turn off internal HDF5 C library error messages.  They are on by default.
    void show_errors(bool on=true);

    /// Ensure all parents to datapath exist as HDF5 groups.
    void ensure_parents(hid_t file_id, const std::string& datapath);

    /** Open or create an HDF5 file
     *
     * @param filename Path to the HDF5 file
     * @param mode File access mode (default: FileMode::rdonly)
     * @return HDF5 file identifier
     * @throws IOError if file cannot be opened or created
     */
    hid_t open(const std::string& filename, unsigned mode = FileMode::rdonly);

    /** Close an HDF5 file
     *
     * @param file_id HDF5 file identifier
     * @throws IOError if file cannot be closed
     */
    void close(hid_t file_id);

    /** Write a dataset from raw bytes
     *
     * @param file_id HDF5 file identifier
     * @param data Pointer to raw data bytes
     * @param shape Dimensions of the dataset
     * @param dtype Data type
     * @param datapath HDF5 path where dataset will be written (e.g., "/group/dataset")
     * @throws IOError on write failure
     */
    void write_dataset(hid_t file_id, const std::byte* data,
                      const std::vector<int64_t>& shape,
                      DataType dtype,
                      const std::string& datapath);

    /** Write a typed dataset (template version)
     *
     * Type is inferred from the vector element type.
     *
     * @param file_id HDF5 file identifier
     * @param data Vector of typed data
     * @param shape Dimensions of the dataset
     * @param datapath HDF5 path where dataset will be written
     * @throws IOError on write failure
     */
    template<typename T>
    void write_dataset(hid_t file_id, const std::vector<T>& data,
                      const std::vector<int64_t>& shape,
                      const std::string& datapath);

    /** Read a dataset into raw bytes
     *
     * @param file_id HDF5 file identifier
     * @param data Output vector to store raw bytes
     * @param shape Output vector to store dimensions
     * @param datapath HDF5 path to dataset (e.g., "/group/dataset")
     * @throws IOError on read failure
     */
    void read_dataset(hid_t file_id, std::vector<std::byte>& data,
                     std::vector<int64_t>& shape,
                     const std::string& datapath);

    /** Read a typed dataset (template version)
     *
     * Type is inferred from the vector element type.
     *
     * @param file_id HDF5 file identifier
     * @param data Output vector to store typed data
     * @param shape Output vector to store dimensions
     * @param datapath HDF5 path to dataset
     * @throws IOError on read failure
     */
    template<typename T>
    void read_dataset(hid_t file_id, std::vector<T>& data,
                     std::vector<int64_t>& shape,
                     const std::string& datapath);

    /** Write metadata as HDF5 attributes
     *
     * Metadata is written as attributes. The datapath can point to either
     * a dataset or a group - the function auto-detects the object type.
     *
     * @param file_id HDF5 file identifier
     * @param metadata Configuration object (Json::Value) containing metadata
     * @param datapath HDF5 path to dataset or group
     * @throws IOError on write failure
     */
    void write_metadata(hid_t file_id, const Configuration& metadata,
                       const std::string& datapath);

    /** Read metadata from HDF5 attributes
     *
     * Reads attributes from the specified datapath. The datapath can point
     * to either a dataset or a group - the function auto-detects.
     *
     * @param file_id HDF5 file identifier
     * @param datapath HDF5 path to dataset or group
     * @return Configuration object containing metadata, empty if no metadata exists
     * @throws IOError on read failure (but not if no attributes exist)
     */
    Configuration read_metadata(hid_t file_id, const std::string& datapath);

    // Explicit instantiation declarations
    extern template void write_dataset<int16_t>(hid_t, const std::vector<int16_t>&,
                                               const std::vector<int64_t>&,
                                               const std::string&);
    extern template void write_dataset<int32_t>(hid_t, const std::vector<int32_t>&,
                                               const std::vector<int64_t>&,
                                               const std::string&);
    extern template void write_dataset<int64_t>(hid_t, const std::vector<int64_t>&,
                                               const std::vector<int64_t>&,
                                               const std::string&);
    extern template void write_dataset<float>(hid_t, const std::vector<float>&,
                                             const std::vector<int64_t>&,
                                             const std::string&);
    extern template void write_dataset<double>(hid_t, const std::vector<double>&,
                                              const std::vector<int64_t>&,
                                              const std::string&);

    extern template void read_dataset<int16_t>(hid_t, std::vector<int16_t>&,
                                              std::vector<int64_t>&,
                                              const std::string&);
    extern template void read_dataset<int32_t>(hid_t, std::vector<int32_t>&,
                                              std::vector<int64_t>&,
                                              const std::string&);
    extern template void read_dataset<int64_t>(hid_t, std::vector<int64_t>&,
                                              std::vector<int64_t>&,
                                              const std::string&);
    extern template void read_dataset<float>(hid_t, std::vector<float>&,
                                            std::vector<int64_t>&,
                                            const std::string&);
    extern template void read_dataset<double>(hid_t, std::vector<double>&,
                                             std::vector<int64_t>&,
                                             const std::string&);

    // Helper to write JSON value as HDF5 attribute
    void write_json_attribute(hid_t loc_id, const std::string& attr_name,
                             const Json::Value& value);

    // Helper to read JSON value from HDF5 attribute
    Json::Value read_json_attribute(hid_t loc_id, const std::string& attr_name);

}  // namespace WireCell::Hio

#endif
