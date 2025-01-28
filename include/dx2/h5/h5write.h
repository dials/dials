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

#endif // H5WRITE_H
