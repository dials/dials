#ifndef H5WRITE_H
#define H5WRITE_H

#include <array>
#include <cstdlib>
#include <hdf5.h>
#include <iostream>
#include <string>
#include <vector>

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

#endif // H5WRITE_H
