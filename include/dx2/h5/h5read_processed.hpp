/**
 * @file h5read_processed.hpp
 * @brief HDF5 dataset access utilities for reflection table loading.
 *
 * This header provides functions for reading typed datasets from HDF5
 * files along with their dimensional metadata. It includes fast,
 * shape-aware readers, recursive traversal of group structures, and
 * scoped utilities for extracting datasets from flat HDF5 groups such
 * as `/dials/processing/group_0`.
 *
 * Utilities in `h5read_processed_utils` focus on immediate dataset
 * access, while recursive tools allow for full introspection of nested
 * groups.
 *
 * Intended for internal use in processing pipelines.
 */

#pragma once

#include "dx2/h5/h5utils.hpp"
#include "dx2/logging.hpp"
#include <cassert>
#include <chrono>
#include <cstring>
#include <experimental/mdspan>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string>
#include <typeindex>
#include <unordered_set>
#include <vector>

using namespace h5utils;

namespace h5read_processed_utils {

#pragma region Internal Structs

/**
 * @brief Container for shaped HDF5 data.
 *
 * Holds the raw flat data vector and shape information for a dataset.
 */
template <typename T> struct H5ArrayData {
  std::vector<T> data;
  std::vector<size_t> shape;
};

/**
 * @brief Internal context for scanning a single HDF5 group.
 *
 * Used to pass state to the `H5Literate2` callback for flat dataset collection.
 */
struct GroupScanContext {
  std::string group_name;            // Name of the group being scanned
  std::vector<std::string> datasets; // Collected dataset paths
};

///  @brief Traverse state context for recursive dataset discovery.
struct TraverseData {
  std::vector<std::string> *datasets;
  std::string current_path;
  std::unordered_set<std::string> *visited_groups;
};

#pragma endregion
#pragma region Callbacks

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

/**
 * @brief Callback for `H5Literate2` that walks both groups and datasets.
 *
 * Recursively visits groups and appends dataset paths to `TraverseData`.
 */
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
    dx2_log::error(fmt::format("Unable to get object info for: {}", full_path));
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
    H5Group group_id(H5Gopen2(loc_id, name, H5P_DEFAULT));
    if (group_id) {
      traverse_hdf5(group_id, full_path, *(traverse_data->datasets),
                    *(traverse_data->visited_groups));
    } else {
      dx2_log::error(fmt::format("Unable to open group: {}", full_path));
    }
  }

  return 0; // Continue iteration
}

#pragma endregion
#pragma region Group Traversal

/**
 * @brief Recursively traverses an HDF5 group and collects dataset paths.
 *
 * @param loc_id HDF5 group ID to start from.
 * @param path Path prefix for entries.
 * @param datasets Output vector to populate.
 * @param visited_groups Prevents revisiting cycles.
 */
void traverse_hdf5(hid_t loc_id, const std::string &path,
                   std::vector<std::string> &datasets,
                   std::unordered_set<std::string> &visited_groups) {
  // std::cout << "Traversing: " << (path.empty() ? "/" : path) << std::endl;

  TraverseData traverse_data = {&datasets, path, &visited_groups};
  H5Literate2(loc_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, group_iterator,
              &traverse_data);
}

#pragma endregion
} // namespace h5read_processed_utils

#pragma region Public API
/**
 * @brief Lists all *immediate* datasets in an HDF5 group (non-recursive).
 *
 * Only returns datasets directly under `group_name`, not nested in subgroups.
 *
 * @param filename Path to the HDF5 file.
 * @param group_name Path to the group (e.g., "/dials/processing/group_0").
 * @return Vector of full dataset paths.
 */
std::vector<std::string> get_datasets_in_group(std::string_view filename,
                                               std::string_view group_name) {
  std::string fname(filename);
  std::string gpath(group_name);
  // Open file
  H5File file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file) {
    throw std::runtime_error("Error: Unable to open file: " + fname);
  }

  // Open group
  H5Group group(H5Gopen2(file, gpath.c_str(), H5P_DEFAULT));
  if (!group) {
    dx2_log::warning(fmt::format("Missing group '{}', skipping.", gpath));
    return {};
  }

  h5read_processed_utils::GroupScanContext context{gpath, {}};

  // Iterate over immediate children
  H5Literate2(group, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
              h5read_processed_utils::scan_group_callback, &context);

  return context.datasets;
}

/**
 * @brief Reads an HDF5 dataset and returns its data and shape.
 *
 * @tparam T Type of data to read (e.g., `double`, `int`).
 * @param filename Path to the HDF5 file.
 * @param dataset_name Full path to the dataset.
 * @return H5ArrayData<T> containing flat data and shape vector.
 */
template <typename T>
h5read_processed_utils::H5ArrayData<T>
read_array_with_shape_from_h5_file(std::string_view filename,
                                   std::string_view dataset_name) {
  auto start_time = std::chrono::high_resolution_clock::now();

  std::string fname(filename);
  std::string dset_name(dataset_name);

  // Open the HDF5 file
  H5File file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file) {
    throw std::runtime_error("Error: Unable to open file: " + fname);
  }

  // Open the dataset
  H5Dataset dataset(H5Dopen(file, dset_name.c_str(), H5P_DEFAULT));
  if (!dataset) {
    throw std::runtime_error("Error: Unable to open dataset: " + dset_name);
  }

  // Get the datatype and check size
  H5Type datatype(H5Dget_type(dataset));
  if (H5Tget_size(datatype) != sizeof(T)) {
    throw std::runtime_error(
        "Error: Dataset type size does not match expected type size.");
  }

  // Get the dataspace and the number of elements
  H5Space dataspace(H5Dget_space(dataset));
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
    throw std::runtime_error("Error: Unable to read dataset: " + dset_name);
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  double elapsed_time =
      std::chrono::duration<double>(end_time - start_time).count();
  dx2_log::debug("READ TIME for {} : {:.4f}s", dset_name, elapsed_time);

  return {std::move(data_out), std::vector<size_t>(dims.begin(), dims.end())};
}

/**
 * @brief Reads an HDF5 dataset into a flat vector (no shape metadata).
 *
 * @tparam T Type of the dataset (e.g., `float`, `int64_t`).
 * @param filename HDF5 file path.
 * @param dataset_name Full path to the dataset.
 * @return Vector of raw values.
 */
template <typename T>
std::vector<T> read_array_from_h5_file(std::string_view filename,
                                       std::string_view dataset_name) {
  return read_array_with_shape_from_h5_file<T>(filename, dataset_name).data;
}

/**
 * @brief Recursively finds all datasets in a group and subgroups.
 *
 * @param filename Path to the HDF5 file.
 * @param group_name Name of the top-level group to search.
 * @return Vector of full dataset paths.
 */
std::vector<std::string>
get_datasets_in_group_recursive(std::string_view filename,
                                std::string_view group_name) {
  std::string fname(filename);
  std::string gpath(group_name);

  std::vector<std::string> datasets;
  std::unordered_set<std::string> visited_groups;

  // Open the HDF5 file
  H5File file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file) {
    dx2_log::error(fmt::format("Unable to open file: {}", fname));
  }

  // Open the group
  // hid_t group = H5Gopen2(file, gpath.c_str(), H5P_DEFAULT);
  H5Group group(H5Gopen2(file, gpath.c_str(), H5P_DEFAULT));
  if (!group) {
    dx2_log::warning(fmt::format("Missing group '{}', skipping.", gpath));
    return {};
  }

  // Start traversal from the group
  h5read_processed_utils::traverse_hdf5(group, gpath, datasets, visited_groups);

  return datasets;
}

inline void read_experiment_metadata(hid_t group_id,
                                     std::vector<uint64_t> &experiment_ids,
                                     std::vector<std::string> &identifiers) {
  if (H5Aexists(group_id, "experiment_ids") > 0) {
    H5Attr attr(H5Aopen(group_id, "experiment_ids", H5P_DEFAULT));
    H5Space space(H5Aget_space(attr));
    hssize_t num_elements = H5Sget_simple_extent_npoints(space);

    experiment_ids.resize(num_elements);
    H5Aread(attr, H5T_NATIVE_ULLONG, experiment_ids.data());
  }

  if (H5Aexists(group_id, "identifiers") > 0) {
    H5Attr attr(H5Aopen(group_id, "identifiers", H5P_DEFAULT));
    H5Type type(H5Aget_type(attr));
    H5Space space(H5Aget_space(attr));
    hssize_t num_elements = H5Sget_simple_extent_npoints(space);

    std::vector<char *> raw_strings(num_elements);
    identifiers.resize(num_elements);
    H5Aread(attr, type, raw_strings.data());

    for (hssize_t i = 0; i < num_elements; ++i) {
      identifiers[i] = std::string(raw_strings[i]);
    }
  }
}

/**
 * @brief Extracts just the leaf name from a full HDF5 dataset path.
 *
 * @param path Full dataset path (e.g., `/a/b/c`).
 * @return The base name (`c`).
 */
std::string get_dataset_name(const std::string &path) {
  size_t pos = path.find_last_of('/');
  if (pos == std::string::npos) {
    return path; // No '/' found, return the whole path
  }
  return path.substr(pos + 1);
}
#pragma endregion