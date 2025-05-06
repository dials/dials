/**
 * @file h5dispatch.hpp
 * @brief Central type registry and dispatch utilities for HDF5-C++ interop.
 *
 * This header provides:
 * - A centralised registry of supported C++/HDF5 type mappings.
 * - Functions to dispatch based on HDF5 dataset types or C++ type_index.
 * - A clean extensible structure to support future type additions.
 */

#pragma once

#include "dx2/h5/h5utils.hpp"
#include <hdf5.h>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <typeindex>
#include <vector>

using namespace h5utils;

namespace h5dispatch {

/**
 * @brief Template wrapper used to tag types in dispatch callbacks.
 *
 * This enables dispatch functions to pass types as compile-time values.
 */
template <typename T> struct type_tag_t {
  using type = T;
};
template <typename T> constexpr type_tag_t<T> type_tag{};

/**
 * @brief Metadata describing a single registered type mapping.
 */
struct H5TypeInfo {
  std::type_index cpp_type; ///< The C++ type_index of the registered type
  hid_t h5_native_type;     ///< The native HDF5 type (e.g. H5T_NATIVE_INT)
  H5T_class_t h5_class;     ///< HDF5 class (integer, float, etc.)
  size_t h5_size;           ///< Size in bytes
  H5T_sign_t h5_sign;       ///< Sign info for integers (signed/unsigned)
};

/**
 * @brief Central registry for all supported type mappings.
 *
 * This is the single source of truth for what C++ types are supported,
 * and how they map to HDF5 native types.
 */
inline const std::vector<H5TypeInfo> &get_supported_types() {
  static const std::vector<H5TypeInfo> registry = {
      {typeid(int), H5T_NATIVE_INT, H5T_INTEGER, sizeof(int), H5T_SGN_2},
      {typeid(int64_t), H5T_NATIVE_LLONG, H5T_INTEGER, sizeof(int64_t),
       H5T_SGN_2},
      {typeid(uint64_t), H5T_NATIVE_ULLONG, H5T_INTEGER, sizeof(uint64_t),
       H5T_SGN_NONE},
      {typeid(double), H5T_NATIVE_DOUBLE, H5T_FLOAT, sizeof(double),
       H5T_SGN_ERROR},
  };
  return registry;
}

/**
 * @brief Returns the native HDF5 type corresponding to a given C++ type.
 *
 * @tparam T The C++ type.
 * @return The HDF5 native type (hid_t).
 * @throws std::runtime_error if the type is not in the registry.
 */
template <typename T> hid_t get_h5_native_type() {
  for (const auto &entry : get_supported_types()) {
    if (entry.cpp_type == typeid(T)) {
      return entry.h5_native_type;
    }
  }
  throw std::runtime_error("Unsupported type in get_h5_native_type().");
}

/**
 * @brief Internal dispatcher used by public-facing dispatchers.
 *
 * Given a matcher function, this tries all registered types and invokes
 * the callback with the matching one.
 *
 * @tparam Callback A callable that takes a type_tag<T>.
 * @tparam MatchFn A predicate function to match an H5TypeInfo entry.
 * @param match_fn The matcher used to find the appropriate type.
 * @param cb The callback to invoke once a match is found.
 *
 * @throws std::runtime_error if no match is found.
 */
template <typename Callback, typename MatchFn>
void dispatch_registered_type(const MatchFn &match_fn, Callback &&cb) {
  for (const auto &entry : get_supported_types()) {
    if (match_fn(entry)) {
      // Dispatch using the appropriate type tag
      if (entry.cpp_type == typeid(int))
        cb(type_tag<int>);
      else if (entry.cpp_type == typeid(int64_t))
        cb(type_tag<int64_t>);
      else if (entry.cpp_type == typeid(uint64_t))
        cb(type_tag<uint64_t>);
      else if (entry.cpp_type == typeid(double))
        cb(type_tag<double>);
      else
        throw std::runtime_error("Type registered but not dispatchable.");
      return;
    }
  }
  throw std::runtime_error("No matching registered type.");
}

/**
 * @brief Dispatch based on the type of an HDF5 dataset.
 *
 * Inspects the HDF5 type of the dataset and invokes the callback
 * with the appropriate type tag.
 *
 * @tparam Callback A callable with templated operator()(type_tag<T>).
 * @param dataset_id The HDF5 dataset to inspect.
 * @param cb The callback to invoke with the matched type tag.
 *
 * @throws std::runtime_error if no supported match is found.
 */
template <typename Callback>
void dispatch_h5_dataset_type(hid_t dataset_id, Callback &&cb) {
  H5Type type_id(H5Dget_type(dataset_id));
  if (!type_id) {
    throw std::runtime_error("Failed to get dataset type.");
  }

  H5T_class_t cls = H5Tget_class(type_id);
  size_t size = H5Tget_size(type_id);
  H5T_sign_t sign = H5Tget_sign(type_id);

  // Matcher lambda for registered types
  auto matcher = [cls, size, sign](const H5TypeInfo &entry) {
    return entry.h5_class == cls && entry.h5_size == size &&
           (cls != H5T_INTEGER || entry.h5_sign == sign);
  };

  try {
    dispatch_registered_type(matcher, std::forward<Callback>(cb));
  } catch (...) {
    std::cerr << "Unsupported HDF5 dataset type:\n"
              << "  class: " << static_cast<int>(cls) << "\n"
              << "  size: " << size << "\n"
              << "  sign: " << static_cast<int>(sign) << "\n";
    throw;
  }
}

/**
 * @brief Dispatch based on a std::type_index (e.g., from a column).
 *
 * Looks up the registered type and invokes the callback with the type tag.
 *
 * @tparam Callback A callable with templated operator()(type_tag<T>).
 * @param t The std::type_index to match.
 * @param cb The callback to invoke with the matched type tag.
 *
 * @throws std::runtime_error if the type is not in the registry.
 */
template <typename Callback>
void dispatch_column_type(const std::type_index &t, Callback &&cb) {
  auto matcher = [t](const H5TypeInfo &entry) { return entry.cpp_type == t; };
  dispatch_registered_type(matcher, std::forward<Callback>(cb));
}

} // namespace h5dispatch