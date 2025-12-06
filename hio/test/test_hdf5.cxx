#include "WireCellUtil/Exceptions.h"

#include <iostream>

#include <hdf5.h>


int main()
{   
    // turn off HDF5-DIAG message globally
    // H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    std::string fn = "g4-rec-0.h5";
    auto hfile_existing = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // seems no way to check if file is already in use, so just try to create it
    hid_t hfile = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (hfile == H5I_INVALID_HID) {
        std::cout << "File " << fn << " exists, opening it instead." << std::endl;
        hfile = H5Fopen(fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    /**
     *  this seems cleaner but not the logic we want
    */
    // if (H5Fis_hdf5(fn.c_str()) > 0) {
    //     std::cout << "File " << fn << " exists, opening it instead." << std::endl;
    //     hfile = H5Fopen(fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // } else {
    //     std::cout << "Creating file " << fn << std::endl;
    //     hfile = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // }

    if (hfile == H5I_INVALID_HID) {
        WireCell::raise<WireCell::IOError>("Failed to create file %s", fn);
    }
    
    // Create a dataset
    int data = 42;
    hid_t dataspace = H5Screate(H5S_SCALAR);
    hid_t dataset = H5Dcreate(hfile, "my_dataset", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(hfile);

    H5Fclose(hfile_existing);

    return 0;
}