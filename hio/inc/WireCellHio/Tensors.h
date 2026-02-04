/** HDF5 serialization for ITensor and ITensorSet
 *
 * Functions to write and read ITensor and ITensorSet objects to/from HDF5 files.
 * Lives in the WireCell::Hio namespace.
 */

#ifndef WIRECELLHIO_TENSORS
#define WIRECELLHIO_TENSORS

#include "WireCellIface/ITensor.h"
#include "WireCellIface/ITensorSet.h"
#include <hdf5.h>
#include <string>

namespace WireCell::Hio {

    /** Write an ITensor to HDF5 file
     *
     * Saves the tensor data as an HDF5 dataset at the specified datapath.
     * The tensor's dtype, shape, and data are preserved. Metadata from
     * ITensor::metadata() is saved as HDF5 attributes on the dataset.
     *
     * @param file_id HDF5 file identifier
     * @param tensor Shared pointer to ITensor to write
     * @param datapath HDF5 path where tensor will be written (e.g., "/group/tensor")
     * @param gzip Compression level 0-9 (0=no compression, 9=max compression)
     * @param chunks User-specified chunk sizes (empty=auto, applies to last dimensions)
     * @throws IOError on write failure or HDF5 errors
     */
    void write_itensor(hid_t file_id, ITensor::pointer tensor, const std::string& datapath,
                      int gzip = 0, const std::vector<int>& chunks = {});

    /** Read an ITensor from HDF5 file
     *
     * Reads tensor data and metadata from HDF5 and constructs a SimpleTensor.
     * The returned tensor will have the same dtype, shape, data, and metadata
     * as was written.
     *
     * @param file_id HDF5 file identifier
     * @param datapath HDF5 path to tensor dataset
     * @return Shared pointer to ITensor (SimpleTensor implementation)
     * @throws IOError on read failure or HDF5 errors
     */
    ITensor::pointer read_itensor(hid_t file_id, const std::string& datapath);

    /** Write an ITensorSet to HDF5 file
     *
     * Saves the tensor set as an HDF5 group at the specified datapath.
     * The set's metadata is saved as attributes on the group. Each tensor
     * in the set is saved as a dataset within the group, named by its index
     * (e.g., "0", "1", "2", etc.) using write_itensor().
     *
     * @param file_id HDF5 file identifier
     * @param tensorset Shared pointer to ITensorSet to write
     * @param datapath HDF5 path where tensor set will be written (e.g., "/group/tensorset")
     * @param gzip Compression level 0-9 (0=no compression, 9=max compression)
     * @param chunks User-specified chunk sizes (empty=auto, applies to last dimensions)
     * @throws IOError on write failure or HDF5 errors
     */
    void write_itensorset(hid_t file_id, ITensorSet::pointer tensorset, const std::string& datapath,
                         int gzip = 0, const std::vector<int>& chunks = {});

    /** Read an ITensorSet from HDF5 file
     *
     * Reads tensor set data and metadata from HDF5 and constructs a SimpleTensorSet.
     * The group's attributes are read as metadata. Each dataset in the group is
     * read as an ITensor using read_itensor(). The returned tensor set will have
     * the same metadata and tensors as was written.
     *
     * @param file_id HDF5 file identifier
     * @param datapath HDF5 path to tensor set group
     * @return Shared pointer to ITensorSet (SimpleTensorSet implementation)
     * @throws IOError on read failure or HDF5 errors
     */
    ITensorSet::pointer read_itensorset(hid_t file_id, const std::string& datapath);

}  // namespace WireCell::Hio

#endif
