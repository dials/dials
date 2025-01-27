#ifndef H5WRITE_H
#define H5WRITE_H

#include <array>
#include <cstdlib>
#include <hdf5.h>
#include <iostream>
#include <string>
#include <vector>

// Function to create a group if it does not exist
hid_t create_or_open_group(hid_t parent, const std::string &group_name) {
  hid_t group = H5Gopen(parent, group_name.c_str(), H5P_DEFAULT);
  if (group < 0) {
    group = H5Gcreate(parent, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (group < 0) {
      std::cerr << "Error: Unable to create group " << group_name << std::endl;
      std::exit(1);
    }
  }
  return group;
}

// Function to write a 3D vector dataset to an HDF5 file
void write_xyzobs_to_h5_file(
    const std::string &filename, const std::string &group_path,
    const std::string &dataset_name,
    const std::vector<std::array<double, 3>> &xyz_data) {
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
      std::cerr << "Error: Unable to create or open file " << filename
                << std::endl;
      std::exit(1);
    }
  }

  // Create or open the hierarchy: /dials/processing/group_0
  hid_t dials_group = create_or_open_group(file, "dials");
  hid_t processing_group = create_or_open_group(dials_group, group_path);
  hid_t group = create_or_open_group(processing_group, "group_0");

  // Create dataspace for 3D vector data
  hsize_t dims[2] = {xyz_data.size(), 3}; // Number of elements x 3D coordinates
  hid_t dataspace = H5Screate_simple(2, dims, NULL);

  // Create dataset
  hid_t dataset = H5Dcreate(group, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    std::cerr << "Error: Unable to create dataset " << dataset_name
              << std::endl;
    std::exit(1);
  }

  // Flatten 3D vector data for HDF5 storage
  std::vector<double> flat_data;
  for (const auto &xyz : xyz_data) {
    flat_data.insert(flat_data.end(), xyz.begin(), xyz.end());
  }

  // Write data to dataset
  herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, flat_data.data());
  if (status < 0) {
    std::cerr << "Error: Unable to write data to dataset " << dataset_name
              << std::endl;
    std::exit(1);
  }

  // Cleanup
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Gclose(group);
  H5Gclose(processing_group);
  H5Gclose(dials_group);
  H5Fclose(file);

  std::cout << "Successfully wrote 3D vector dataset " << dataset_name << " to "
            << filename << std::endl;
}

#endif // H5WRITE_H
