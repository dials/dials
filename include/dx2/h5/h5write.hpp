/**
 * @file h5write.hpp
 * @brief Utilities for writing raw and structured data to HDF5 files.
 *
 * This header provides two levels of HDF5 writing utilities:
 *
 * - **Raw Writer:** Directly writes typed raw memory buffers into
 *   datasets. Suitable for low-level, performance-sensitive operations
 *   like reflection table saving.
 *
 * - **High-Level Writer:** Dynamically infers shape from nested
 *   standard containers (e.g., `std::vector<std::vector<double>>`),
 *   flattens them, and writes them to HDF5. Useful for quickly
 *   serializing in-memory C++ containers to disk.
 *
 * Intended for internal use in DX2 processing and general scientific
 * data workflows.
 */
#pragma once

#include <array>
#include <cstdlib>
#include <hdf5.h>
#include <iostream>
#include <string>
#include <vector>

#pragma region Raw writer
/**
 * @brief Recursively traverses or creates groups in an HDF5 file based
 * on the given path.
 *
 * This function takes a parent group identifier and a path string, and
 * recursively traverses or creates the groups specified in the path. If
 * a group in the path does not exist, it is created.
 *
 * @param parent The identifier of the parent group in the HDF5 file.
 * @param path The path of groups to traverse or create, specified as a
 * string with '/' as the delimiter.
 * @return The identifier of the final group in the path.
 * @throws std::runtime_error If a group cannot be created or opened.
 */
hid_t traverse_or_create_groups(hid_t parent, const std::string &path) {
  // Strip leading '/' characters, if any, to prevent empty group names
  size_t start_pos = path.find_first_not_of('/');
  if (start_pos == std::string::npos) {
    return parent; // Return parent if the path is entirely '/'
  }
  std::string cleaned_path = path.substr(start_pos);

  /*
   * This is the base case for recursion. When the path is empty, we
   * have reached the final group in the path and we return the parent
   * group.
   */
  if (cleaned_path.empty()) {
    return parent;
  }

  // Split the path into the current group name and the remaining path
  size_t pos = cleaned_path.find('/');
  std::string group_name =
      (pos == std::string::npos) ? cleaned_path : cleaned_path.substr(0, pos);
  std::string remaining_path =
      (pos == std::string::npos) ? "" : cleaned_path.substr(pos + 1);

  // Attempt to open the group. If it does not exist, create it.
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL); // Suppress errors to stdout when
  // trying to open a file/group that may not exist.
  hid_t next_group = H5Gopen(parent, group_name.c_str(), H5P_DEFAULT);
  if (next_group < 0) {
    next_group = H5Gcreate(parent, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
    if (next_group < 0) {
      std::runtime_error("Error: Unable to create or open group: " +
                         group_name);
    }
  }

  // Recurse to the next group in the hierarchy
  hid_t final_group = traverse_or_create_groups(next_group, remaining_path);

  // Close the current group to avoid resource leaks, except for the final group
  if (next_group != final_group) {
    H5Gclose(next_group);
  }

  return final_group;
}

/**
 * @brief Writes raw data to an HDF5 file.
 *
 * This function writes a dataset to an HDF5 file. The dataset's shape
 * is specified by the user.
 *
 * @param filename The path to the HDF5 file.
 * @param dataset_path The full path to the dataset, including group
 * hierarchies.
 * @param data_ptr Pointer to the raw data to write to the dataset.
 * @param shape The shape of the dataset as a vector of dimensions.
 * @throws std::runtime_error If the dataset cannot be created or data
 * cannot be written.
 */
template <typename T>
void write_raw_data_to_h5_file(const std::string &filename,
                               const std::string &dataset_path,
                               const T *data_ptr,
                               const std::vector<hsize_t> &shape) {
  // Open or create the file
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
      throw std::runtime_error("Failed to create or open file: " + filename);
    }
  }

  try {
    // Split group and dataset name
    size_t slash_pos = dataset_path.find_last_of('/');
    if (slash_pos == std::string::npos) {
      throw std::runtime_error("Invalid dataset path: " + dataset_path);
    }
    std::string group_path = dataset_path.substr(0, slash_pos);
    std::string dataset_name = dataset_path.substr(slash_pos + 1);

    hid_t group = traverse_or_create_groups(file, group_path);

    // Create dataspace
    hid_t dataspace = H5Screate_simple(shape.size(), shape.data(), nullptr);

    // Determine native HDF5 type
    hid_t h5_type;
    if constexpr (std::is_same_v<T, double>)
      h5_type = H5T_NATIVE_DOUBLE;
    else if constexpr (std::is_same_v<T, int>)
      h5_type = H5T_NATIVE_INT;
    else if constexpr (std::is_same_v<T, int64_t>)
      h5_type = H5T_NATIVE_LLONG;
    else if constexpr (std::is_same_v<T, uint64_t>)
      h5_type = H5T_NATIVE_ULLONG;
    else
      throw std::runtime_error("Unsupported type for HDF5 writing");

    // Create or overwrite dataset
    hid_t dset = H5Dcreate2(group, dataset_name.c_str(), h5_type, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
      dset = H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT);
      if (dset < 0) {
        H5Sclose(dataspace);
        throw std::runtime_error("Failed to create or open dataset: " +
                                 dataset_name);
      }
    }

    herr_t status =
        H5Dwrite(dset, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_ptr);
    if (status < 0) {
      throw std::runtime_error("Failed to write dataset: " + dataset_path);
    }

    H5Dclose(dset);
    H5Sclose(dataspace);
    H5Gclose(group);
  } catch (...) {
    H5Fclose(file);
    throw;
  }

  H5Fclose(file);
}
#pragma endregion
#pragma region High-level writer
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The high-level writer is useful for writing nested containers. While
 * not utilised within the DX2 backend itself, it is a useful utility
 * for writing data to HDF5 files, thus it is included here.
 */

/**
 * @brief Deduce the shape of a nested container.
 *
 * This helper function recursively determines the shape of nested
 * containers. This is used to determine the shape of the dataset to be
 * created in an HDF5 file.
 *
 * @tparam Container The type of the container.
 * @param container The container whose shape is to be determined.
 * @return A vector of dimensions representing the shape of the
 * container.
 */
template <typename Container>
std::vector<hsize_t> deduce_shape(const Container &container) {
  if (container.empty()) {
    return {0};
  }
  if constexpr (std::is_arithmetic_v<typename Container::value_type>) {
    // Base case: container holds arithmetic types (e.g., double, int)
    return {container.size()};
  } else {
    // Recursive case: container holds other containers

    // Check that all inner containers have the same size
    size_t inner_size = container.begin()->size();
    for (const auto &sub_container : container) {
      if (sub_container.size() != inner_size) {
        throw std::runtime_error("Cannot deduce shape: inner containers have "
                                 "different sizes.");
      }
    }

    auto sub_shape = deduce_shape(*container.begin());
    sub_shape.insert(sub_shape.begin(), container.size());
    return sub_shape;
  }
}

/**
 * @brief Flatten nested containers into a 1D vector.
 *
 * This helper function recursively flattens nested containers into a 1D
 * vector for writing to HDF5. If the input container is already 1D, it
 * simply returns it.
 *
 * @tparam Container The type of the container.
 * @param container The container to flatten.
 * @return A flat vector containing all elements of the input container
 * in a linear order.
 */
template <typename Container> auto flatten(const Container &container) {
  // Determine the type of the elements in the container
  using ValueType = typename Container::value_type;

  // Base case: If the container holds arithmetic types (e.g., int, double),
  // it is already 1D, so we return a copy of the container as a std::vector.
  if constexpr (std::is_arithmetic_v<ValueType>) {
    return std::vector<ValueType>(container.begin(), container.end());
  } else {
    // Recursive case: The container holds nested containers, so we need to
    // flatten them.

    // Determine the type of elements in the flattened result.
    // This is deduced by recursively calling flatten on the first
    // sub-container.
    using InnerType =
        typename decltype(flatten(*container.begin()))::value_type;

    // Create a vector to store the flattened data
    std::vector<InnerType> flat_data;

    // Iterate over the outer container
    for (const auto &sub_container : container) {
      // Recursively flatten each sub-container
      auto sub_flat = flatten(sub_container);

      // Append the flattened sub-container to the result
      flat_data.insert(flat_data.end(), sub_flat.begin(), sub_flat.end());
    }

    // Return the fully flattened data
    return flat_data;
  }
}

/**
 * @brief Writes multidimensional data to an HDF5 file.
 *
 * This function writes a dataset to an HDF5 file. The dataset's shape
 * is determined dynamically based on the input container.
 *
 * @tparam Container The type of the container holding the data.
 * @param filename The path to the HDF5 file.
 * @param dataset_path The full path to the dataset, including group
 * hierarchies.
 * @param data The data to write to the dataset.
 * @throws std::runtime_error If the dataset cannot be created or data
 * cannot be written.
 */
template <typename Container>
void write_data_to_h5_file(const std::string &filename,
                           const std::string &dataset_path,
                           const Container &data) {
  // Open or create the HDF5 file
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL); // Suppress errors to stdout when
  // trying to open a file/group that may not exist.
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
      throw std::runtime_error("Error: Unable to create or open file: " +
                               filename);
    }
  }

  try {
    // Separate the dataset path into group path and dataset name
    size_t last_slash_pos = dataset_path.find_last_of('/');
    if (last_slash_pos == std::string::npos) {
      throw std::runtime_error("Error: Invalid dataset path, no '/' found: " +
                               dataset_path);
    }

    std::string group_path = dataset_path.substr(0, last_slash_pos);
    std::string dataset_name = dataset_path.substr(last_slash_pos + 1);

    // Traverse or create the groups leading to the dataset
    hid_t group = traverse_or_create_groups(file, group_path);

    // Deduce the shape of the data
    std::vector<hsize_t> shape = deduce_shape(data);

    // Flatten the data into a 1D vector
    auto flat_data = flatten(data);

    // Check if dataset exists
    hid_t dataset = H5Dopen(group, dataset_name.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
      // Dataset does not exist, create it
      hid_t dataspace = H5Screate_simple(shape.size(), shape.data(), NULL);
      if (dataspace < 0) {
        throw std::runtime_error(
            "Error: Unable to create dataspace for dataset: " + dataset_name);
      }

      dataset = H5Dcreate(group, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                          dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (dataset < 0) {
        H5Sclose(dataspace);
        throw std::runtime_error("Error: Unable to create dataset: " +
                                 dataset_name);
      }

      H5Sclose(dataspace);
    } else {
      // Dataset exists, check if the shape matches
      hid_t existing_space = H5Dget_space(dataset);
      int ndims = H5Sget_simple_extent_ndims(existing_space);
      std::vector<hsize_t> existing_dims(ndims);
      H5Sget_simple_extent_dims(existing_space, existing_dims.data(), NULL);
      H5Sclose(existing_space);

      if (existing_dims != shape) {
        H5Dclose(dataset);
        throw std::runtime_error(
            "Error: Dataset shape mismatch. Cannot overwrite dataset: " +
            dataset_name);
      }

      // Dataset exists and has the correct shape, proceed to overwrite
    }

    // Write the data to the dataset
    herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, flat_data.data());
    if (status < 0) {
      H5Dclose(dataset);
      throw std::runtime_error("Error: Unable to write data to dataset: " +
                               dataset_name);
    }

    // Cleanup resources
    H5Dclose(dataset);
    H5Gclose(group);
  } catch (...) {
    H5Fclose(file);
    throw; // Re-throw the exception to propagate it upwards
  }

  // Close the file
  H5Fclose(file);
}
#pragma endregion
#pragma region Attribute writer
inline void
write_experiment_metadata(hid_t group_id,
                          const std::vector<uint64_t> &experiment_ids,
                          const std::vector<std::string> &identifiers) {
  // Write experiment_ids
  {
    hsize_t dims = experiment_ids.size();
    hid_t space = H5Screate_simple(1, &dims, nullptr);
    hid_t attr = H5Acreate2(group_id, "experiment_ids", H5T_NATIVE_ULLONG,
                            space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_ULLONG, experiment_ids.data());
    H5Aclose(attr);
    H5Sclose(space);
  }

  // Write identifiers
  {
    hsize_t dims = identifiers.size();
    hid_t space = H5Screate_simple(1, &dims, nullptr);

    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, H5T_VARIABLE);
    H5Tset_cset(str_type, H5T_CSET_UTF8);
    H5Tset_strpad(str_type, H5T_STR_NULLTERM);

    std::vector<const char *> c_strs;
    for (const auto &s : identifiers) {
      c_strs.push_back(s.c_str());
    }

    hid_t attr = H5Acreate2(group_id, "identifiers", str_type, space,
                            H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, str_type, c_strs.data());

    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(str_type);
  }
}
#pragma endregion