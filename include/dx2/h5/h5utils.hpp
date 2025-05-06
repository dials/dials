/**
 * @file h5utils.hpp
 * @brief Provides RAII wrappers for HDF5 resources to ensure proper
 * cleanup and type safety.
 *
 * This header defines a template-based RAII utility, `H5Cleanup`, for
 * managing HDF5 resource lifetimes. It ensures that HDF5 objects are
 * properly closed when they go out of scope, preventing resource leaks.
 * Additionally, it provides type aliases for common HDF5 object types
 * such as files, groups, datasets, attributes, dataspaces, and
 * datatypes.
 *
 * RAII - Resource Acquisition Is Initialization
 * Where resource management is tied to object lifetime.
 */

#pragma once

#include <hdf5.h>
#include <utility>

namespace h5utils {

/// RAII wrapper for HDF5 resources to ensure proper cleanup
template <herr_t(D)(hid_t)> struct H5Cleanup {
  H5Cleanup() : id(-1) {}
  explicit H5Cleanup(hid_t id) : id(id) {}
  H5Cleanup(const H5Cleanup &) = delete;
  H5Cleanup(H5Cleanup &&other) noexcept : id(other.id) { other.id = -1; }
  H5Cleanup &operator=(H5Cleanup &&other) noexcept {
    std::swap(id, other.id);
    return *this;
  }
  ~H5Cleanup() {
    if (id >= 0) {
      D(id);
    }
  }

  operator hid_t() const { return id; }
  explicit operator bool() const { return id >= 0; }
  hid_t id;
};

// Common aliases for specific HDF5 object types
using H5File = H5Cleanup<H5Fclose>;
using H5Group = H5Cleanup<H5Gclose>;
using H5Dataset = H5Cleanup<H5Dclose>;
using H5Attr = H5Cleanup<H5Aclose>;
using H5Space = H5Cleanup<H5Sclose>;
using H5Type = H5Cleanup<H5Tclose>;

} // namespace h5utils