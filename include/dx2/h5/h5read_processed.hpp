#pragma once

#include <cassert>
#include <chrono>
#include <cstring>
#include <experimental/mdspan>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <string>
#include <typeindex>
#include <unordered_set>
#include <vector>

namespace h5read_processed_utils {
/**
 * @brief Container for shaped HDF5 data.
 */
template <typename T> struct H5ArrayData {
  std::vector<T> data;
  std::vector<size_t> shape;
};

/**
 * @brief Internal context for scanning a single HDF5 group.
 *
 * Used to pass state to the `H5Literate2` callback function.
 */
struct GroupScanContext {
  std::string group_name;            /// Name of the group being scanned
  std::vector<std::string> datasets; /// Collected dataset paths
};

/**
 * @brief Callback for `H5Literate2` that collects only datasets.
 *
 * Appends fully-qualified dataset names to the context, skipping subgroups.
 *
 * @param loc_id Current group/location ID.
 * @param name Name of the object (dataset or group).
 * @param info Unused link info.
 * @param op_data Pointer to a `GroupScanContext`.
 * @return 0 on success, non-zero to stop iteration.
 */
herr_t scan_group_callback(hid_t loc_id, const char *name, const H5L_info2_t *,
                           void *op_data) {
  auto *context = static_cast<GroupScanContext *>(op_data);

  // Get object info
  H5O_info2_t obj_info;
  if (H5Oget_info_by_name3(loc_id, name, &obj_info, H5O_INFO_BASIC,
                           H5P_DEFAULT) < 0) {
    return 0; // Skip unreadable entries
  }

  if (obj_info.type == H5O_TYPE_DATASET) {
    context->datasets.push_back(context->group_name + "/" + std::string(name));
  }

  return 0;
}

// Forward declaration
void traverse_hdf5(hid_t loc_id, const std::string &path,
                   std::vector<std::string> &datasets,
                   std::unordered_set<std::string> &visited_groups);

struct TraverseData {
  std::vector<std::string> *datasets;
  std::string current_path;
  std::unordered_set<std::string> *visited_groups;
};

// Callback function for iterating through HDF5 objects
herr_t group_iterator(hid_t loc_id, const char *name, const H5L_info2_t *info,
                      void *op_data) {
  TraverseData *traverse_data = static_cast<TraverseData *>(op_data);
  std::string full_path =
      traverse_data->current_path.empty()
          ? "/" + std::string(name)
          : traverse_data->current_path + "/" + std::string(name);

  // std::cout << "Visiting: " << full_path << std::endl;

  // Ensure this path is not visited twice
  if (traverse_data->visited_groups->find(full_path) !=
      traverse_data->visited_groups->end()) {
    // std::cout << "Skipping already visited group: " << full_path <<
    // std::endl;
    return 0; // Continue iteration
  }

  // Get object info
  H5O_info2_t obj_info;
  if (H5Oget_info_by_name3(loc_id, name, &obj_info, H5O_INFO_BASIC,
                           H5P_DEFAULT) < 0) {
    std::cerr << "Error: Unable to get object info for: " << full_path
              << std::endl;
    return -1;
  }

  if (obj_info.type == H5O_TYPE_DATASET) {
    // std::cout << "Dataset found: " << full_path << std::endl;
    traverse_data->datasets->push_back(full_path);
  } else if (obj_info.type == H5O_TYPE_GROUP) {
    // std::cout << "Entering group: " << full_path << std::endl;

    // Mark the group as visited
    traverse_data->visited_groups->insert(full_path);

    // Open the group to recurse into it
    hid_t group_id = H5Gopen2(loc_id, name, H5P_DEFAULT);
    if (group_id >= 0) {
      traverse_hdf5(group_id, full_path, *(traverse_data->datasets),
                    *(traverse_data->visited_groups));
      H5Gclose(group_id);
    } else {
      std::cerr << "Error: Unable to open group: " << full_path << std::endl;
    }
  }

  return 0; // Continue iteration
}

void traverse_hdf5(hid_t loc_id, const std::string &path,
                   std::vector<std::string> &datasets,
                   std::unordered_set<std::string> &visited_groups) {
  // std::cout << "Traversing: " << (path.empty() ? "/" : path) << std::endl;

  TraverseData traverse_data = {&datasets, path, &visited_groups};
  H5Literate2(loc_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, group_iterator,
              &traverse_data);
}

} // namespace h5read_processed_utils

/**
 * @brief Lists all *immediate* datasets in an HDF5 group (non-recursive).
 *
 * Does not include datasets inside subgroups.
 *
 * @param filename Path to the HDF5 file.
 * @param group_name Path to the group to scan (e.g.,
 * "/dials/processing/group_0").
 * @return Vector of full dataset paths.
 */
std::vector<std::string> get_datasets_in_group(const std::string &filename,
                                               const std::string &group_name) {
  // Open file
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("Error: Unable to open file: " + filename);
  }

  // Open group
  hid_t group = H5Gopen2(file, group_name.c_str(), H5P_DEFAULT);
  if (group < 0) {
    H5Fclose(file);
    std::cerr << "Warning: Missing group " << group_name << ", skipping.\n";
    return {};
  }

  h5read_processed_utils::GroupScanContext context{group_name, {}};

  // Iterate over immediate children
  H5Literate2(group, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
              h5read_processed_utils::scan_group_callback, &context);

  H5Gclose(group);
  H5Fclose(file);

  return context.datasets;
}

/**
 * @brief Reads a dataset from an HDF5 file into a vector with shape.
 *
 * This function reads a dataset from an HDF5 file and returns the data
 * along with its shape. It handles the opening and closing of the
 * HDF5 file and dataset, and checks for errors during the process.
 *
 * @tparam T The type of data to read (e.g., int, double).
 * @param filename The path to the HDF5 file.
 * @param dataset_name The name of the dataset to read.
 * @return A H5ArrayData<T> containing the data and its shape.
 * @throws std::runtime_error If the file, dataset, or datatype cannot be opened
 * or read.
 */
template <typename T>
h5read_processed_utils::H5ArrayData<T>
read_array_with_shape_from_h5_file(const std::string &filename,
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
      int ndims = H5Sget_simple_extent_ndims(dataspace);
      if (ndims <= 0) {
        throw std::runtime_error("Error: Dataset has invalid dimensionality.");
      }

      // Get the dimensions of the dataset
      std::vector<hsize_t> dims(ndims);
      H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);

      // Calculate the number of elements
      size_t num_elements = 1;
      for (auto d : dims)
        num_elements *= d;

      // Allocate a vector to hold the data
      std::vector<T> data_out(num_elements);

      // Read the data into the vector
      herr_t status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              data_out.data());
      if (status < 0) {
        throw std::runtime_error("Error: Unable to read dataset: " +
                                 dataset_name);
      }

      H5Sclose(dataspace);
      H5Tclose(datatype);
      H5Dclose(dataset);
      H5Fclose(file);

      auto end_time = std::chrono::high_resolution_clock::now();
      double elapsed_time =
          std::chrono::duration<double>(end_time - start_time).count();
      std::cout << "READ TIME for " << dataset_name << " : " << elapsed_time
                << "s" << std::endl;

      return {std::move(data_out),
              std::vector<size_t>(dims.begin(), dims.end())};

    } catch (...) {
      H5Dclose(dataset);
      throw;
    }

  } catch (...) {
    H5Fclose(file);
    throw;
  }
}

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
  return read_array_with_shape_from_h5_file<T>(filename, dataset_name).data;
}

/**
 * @brief Discovers all datasets in a group in an HDF5 file.
 *
 * @param filename The path to the HDF5 file.
 * @param group_name The name of the group to search.
 * @return A vector containing the dataset paths.
 */
std::vector<std::string>
get_datasets_in_group_recursive(const std::string &filename,
                                const std::string &group_name) {
  std::vector<std::string> datasets;
  std::unordered_set<std::string> visited_groups;

  // Open the HDF5 file
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("Error: Unable to open file: " + filename);
  }

  // Open the group
  hid_t group = H5Gopen2(file, group_name.c_str(), H5P_DEFAULT);
  if (group < 0) {
    H5Fclose(file);
    std::cerr << "Warning: Missing group " << group_name << ", skipping.\n";
    return {};
  }

  // Start traversal from the group
  h5read_processed_utils::traverse_hdf5(group, group_name, datasets,
                                        visited_groups);

  // Close the group and file
  H5Gclose(group);
  H5Fclose(file);

  return datasets;
}

/**
 * @brief Extracts the dataset name from a full path.
 *
 * @param path The full path to the dataset.
 * @return The name of the dataset.
 */
std::string get_dataset_name(const std::string &path) {
  size_t pos = path.find_last_of('/');
  if (pos == std::string::npos) {
    return path; // No '/' found, return the whole path
  }
  return path.substr(pos + 1);
}
