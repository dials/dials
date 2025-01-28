#ifndef DX2_H5_H5READ_PROCESSED_HPP
#define DX2_H5_H5READ_PROCESSED_HPP

#include <cassert>
#include <chrono>
#include <cstring>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Reads a dataset from an HDF5 file into an std::vector.
 *
 * @tparam T The type of data to read (e.g., int, double).
 * @param filename The path to the HDF5 file.
 * @param dataset_name The name of the dataset to read.
 * @return A std::vector containing the data from the dataset.
 * @throws std::runtime_error If the file, dataset, or datatype cannot be opened
 * or read.
 */
template <typename T>
std::vector<T> read_array_from_h5_file(const std::string &filename,
                                       const std::string &dataset_name) {
  // Start measuring time
  auto start_time = std::chrono::high_resolution_clock::now();

  // Open the HDF5 file
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("Error: Unable to open file: " + filename);
  }

  try {
    // Open the dataset
    hid_t dataset = H5Dopen(file, dataset_name.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
      throw std::runtime_error("Error: Unable to open dataset: " +
                               dataset_name);
    }

    try {
      // Get the datatype and check size
      hid_t datatype = H5Dget_type(dataset);
      size_t datatype_size = H5Tget_size(datatype);
      if (datatype_size != sizeof(T)) {
        throw std::runtime_error(
            "Error: Dataset type size does not match expected type size.");
      }

      // Get the dataspace and the number of elements
      hid_t dataspace = H5Dget_space(dataset);
      size_t num_elements = H5Sget_simple_extent_npoints(dataspace);

      // Allocate a vector to hold the data
      std::vector<T> data_out(num_elements);

      // Read the data into the vector
      herr_t status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              data_out.data());
      if (status < 0) {
        throw std::runtime_error("Error: Unable to read dataset: " +
                                 dataset_name);
      }

      // Close the dataset and return the data
      H5Dclose(dataset);
      H5Fclose(file);

      // Log timing
      auto end_time = std::chrono::high_resolution_clock::now();
      double elapsed_time =
          std::chrono::duration<double>(end_time - start_time).count();
      std::cout << "READ TIME for " << dataset_name << " : " << elapsed_time
                << "s" << std::endl;

      return data_out;

    } catch (...) {
      H5Dclose(dataset); // Ensure dataset is closed in case of failure
      throw;
    }
  } catch (...) {
    H5Fclose(file); // Ensure file is closed in case of failure
    throw;
  }
}

#endif // DX2_H5_H5READ_PROCESSED_HPP