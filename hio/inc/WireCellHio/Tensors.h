/** HDF5 serialization for ITensor
 *
 * Functions to write and read ITensor objects to/from HDF5 files.
 * Lives in the WireCell::Hio namespace.
 */

#ifndef WIRECELLHIO_TENSORS
#define WIRECELLHIO_TENSORS

#include "WireCellIface/ITensor.h"
#include <hdf5.h>
#include <string>

namespace WireCell::Hio {

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

    /** Write an ITensor to HDF5 file
     *
     * Saves the tensor data as an HDF5 dataset at the specified datapath.
     * The tensor's dtype, shape, and data are preserved. Metadata from
     * ITensor::metadata() is saved as HDF5 attributes on the dataset.
     *
     * @param file_id HDF5 file identifier
     * @param tensor Shared pointer to ITensor to write
     * @param datapath HDF5 path where tensor will be written (e.g., "/group/tensor")
     * @throws IOError on write failure or HDF5 errors
     */
    void write_itensor(hid_t file_id, ITensor::pointer tensor, const std::string& datapath);

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

}  // namespace WireCell::Hio

#endif
