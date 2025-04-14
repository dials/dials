#pragma once

#include "dx2/h5/h5read_processed.hpp"
#include "dx2/h5/h5write.hpp"
#include <experimental/mdspan>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#pragma region Utils
/*
 * Type aliases for centralised control of storage types.
 */

// Internal storage type for the data
template <typename T> using storage_type = std::vector<T>;
// Type alias for mdspan (currently using backport std::experimental::mdspan)
template <typename T>
using mdspan_type =
    std::experimental::mdspan<T, std::experimental::dextents<size_t, 2>>;

namespace reflection_table_type_utils {
template <typename T> struct type_tag_t {
  using type = T;
};
template <typename T> constexpr type_tag_t<T> type_tag{};
} // namespace reflection_table_type_utils

/**
 * @brief Dispatches based on the type of an HDF5 dataset.
 *
 * This function inspects the type of the dataset identified by
 * `dataset_id`, and then calls the provided callback as a templated
 * function with the appropriate C++ type.
 *
 * @tparam Callback A callable with a templated `operator()<T>()`.
 * @param dataset_id An open HDF5 dataset identifier (hid_t).
 * @param cb The callback to invoke with the deduced template type.
 */
template <typename Callback>
void dispatch_h5_dataset_type(hid_t dataset_id, Callback &&cb) {
  using namespace reflection_table_type_utils;

  // Get the dataset's datatype
  hid_t type_id = H5Dget_type(dataset_id);
  // Get the class of the datatype (e.g., integer, float)
  H5T_class_t cls = H5Tget_class(type_id);
  // Get the size of the datatype
  size_t size = H5Tget_size(type_id);
  // Get the order of the datatype (e.g., little-endian, big-endian)
  H5T_order_t order = H5Tget_order(type_id);
  // Get the sign of the datatype (e.g., signed, unsigned)
  H5T_sign_t sign = H5Tget_sign(type_id);

  // Match to supported types
  if (cls == H5T_FLOAT && size == sizeof(double)) {
    cb(type_tag<double>);
  } else if (cls == H5T_INTEGER && size == sizeof(int) && sign == H5T_SGN_2) {
    cb(type_tag<int>);
  } else if (cls == H5T_INTEGER && size == sizeof(int64_t) &&
             sign == H5T_SGN_2) {
    cb(type_tag<int64_t>);
  } else if (cls == H5T_INTEGER && size == sizeof(uint64_t) &&
             sign == H5T_SGN_NONE) {
    cb(type_tag<uint64_t>);
  } else {
    // Print full diagnostic
    std::cerr << "Unsupported dataset type:\n";
    std::cerr << "  class: " << static_cast<int>(cls) << "\n";
    std::cerr << "  size: " << size << "\n";
    std::cerr << "  order: " << static_cast<int>(order) << "\n";
    std::cerr << "  sign: " << static_cast<int>(sign) << "\n";
    throw std::runtime_error("Unsupported HDF5 dataset type.");
  }

  H5Tclose(type_id);
}

/**
 * @brief Typeless container for HDF5 data.
 *
 * This is a base class for HDF5 data types. It provides a common interface
 * for different data types, allowing them to be handled polymorphically.
 * Derived classes should implement the virtual methods to provide
 * specific functionality for each data type.
 */
struct ColumnBase {
  virtual ~ColumnBase() = default;
  virtual std::string get_name() const = 0;
  virtual std::vector<size_t> get_shape() const = 0;
  virtual std::type_index get_type() const = 0;
  virtual std::unique_ptr<ColumnBase>
  clone_filtered(const std::vector<size_t> &selected_rows) const = 0;
};

/**
 * @brief Specialized container for HDF5 data of a specific type.
 */
template <typename T> struct TypedColumn : public ColumnBase {
  TypedColumn(const std::string &name, const std::vector<size_t> &shape)
      : name(name), shape(shape), data(compute_size(shape)) {
    setup_mdspan();
  }

  TypedColumn(const std::string &name, const std::vector<size_t> &shape,
              const std::vector<T> &data)
      : name(name), shape(shape), data(data) {
    setup_mdspan();
  }

  TypedColumn(const std::string &name, size_t rows, size_t cols)
      : name(name), shape{rows, cols}, data(rows * cols) {
    setup_mdspan();
  }

  std::vector<size_t> get_shape() const override { return shape; }
  std::string get_name() const override { return name; }
  std::type_index get_type() const override {
    return std::type_index(typeid(T));
  }

  const mdspan_type<T> &span() const { return mdspan_data; }

  std::unique_ptr<ColumnBase>
  clone_filtered(const std::vector<size_t> &selected_rows) const override {
    size_t new_rows = selected_rows.size();
    size_t cols = shape.size() > 1 ? shape[1] : 1;

    std::vector<T> new_data;
    new_data.reserve(new_rows * cols);

    for (size_t row : selected_rows) {
      for (size_t col = 0; col < cols; ++col) {
        new_data.push_back(mdspan_data(row, col));
      }
    }

    return std::make_unique<TypedColumn<T>>(
        name, std::vector<size_t>{new_rows, cols}, new_data);
  }

  std::string name;           // Name of the dataset
  std::vector<size_t> shape;  // Shape of the dataset
  storage_type<T> data;       // Storage for the dataset
  mdspan_type<T> mdspan_data; // mdspan view of the data

private:
  void setup_mdspan() {
    if (shape.size() == 1) {
      mdspan_data = mdspan_type<T>(data.data(), shape[0], 1);
    } else if (shape.size() == 2) {
      mdspan_data = mdspan_type<T>(data.data(), shape[0], shape[1]);
    } else {
      throw std::runtime_error("TypedColumn: unsupported rank " +
                               std::to_string(shape.size()));
    }
  }

  static size_t compute_size(const std::vector<size_t> &shape) {
    if (shape.empty())
      return 0;
    size_t s = 1;
    for (auto dim : shape)
      s *= dim;
    return s;
  }
};
#pragma endregion

#pragma region ReflectionTable
class ReflectionTable {
private:
  /*
   * Stores data in a typeless, polymorphic container
   */
  std::vector<std::unique_ptr<ColumnBase>> data;
  std::string h5_filepath;
  const std::string DEFAULT_REFL_GROUP = "/dials/processing/group_0";

  /// Get the number of rows in the first column, assuming all columns
  /// have the same number of rows
  size_t get_row_count() const {
    if (data.empty())
      return 0;
    return data.front()->get_shape()[0];
  }

public:
  ReflectionTable() = default;

  ReflectionTable(const std::string &h5_filepath) : h5_filepath(h5_filepath) {
    // Discover all datasets in the default reflection group
    std::vector<std::string> datasets =
        get_datasets_in_group(h5_filepath, DEFAULT_REFL_GROUP);

    // Open the HDF5 file once for all dataset operations
    hid_t file = H5Fopen(h5_filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      throw std::runtime_error("Could not open file: " + h5_filepath);
    }

    // Loop over every dataset path in the group
    for (const auto &dataset : datasets) {
      std::string dataset_name = get_dataset_name(dataset);

      // Open the specific dataset (within the file opened above)
      hid_t dataset_id = H5Dopen(file, dataset.c_str(), H5P_DEFAULT);
      if (dataset_id < 0) {
        std::cerr << "Could not open dataset: " << dataset << "\n";
        continue;
      }

      try {
        // Dispatch based on the dataset's actual HDF5 type (float64, int32,
        // etc.)
        dispatch_h5_dataset_type(dataset_id, [&](auto tag) {
          // Deduce the type
          using T = typename decltype(tag)::type;

          // Read the dataset as type T (returns raw data and shape)
          auto result =
              read_array_with_shape_from_h5_file<T>(h5_filepath, dataset);

          // Create a TypedColumn<T> and store it in the typeless container
          data.push_back(std::make_unique<TypedColumn<T>>(
              dataset_name, result.shape, result.data));

          std::cout << "Loaded column: " << dataset_name << " with type "
                    << typeid(T).name() << "\n";
        });
      } catch (const std::exception &e) {
        // If dispatch or reading fails, skip this dataset
        std::cerr << "Skipping dataset " << dataset << ": " << e.what() << "\n";
      }

      H5Dclose(dataset_id);
    }

    H5Fclose(file);
  }

  /**
   * @brief Get a list of all column names in the table.
   */
  std::vector<std::string> get_column_names() const {
    std::vector<std::string> names;
    for (const auto &col : data) {
      names.push_back(col->get_name());
    }
    return names;
  }

  /**
   * @brief Get a typed column by name. Returns nullptr if name or type
   * mismatch.
   *
   * @tparam T The expected type of the column.
   * @param name The dataset/column name.
   * @return A pointer to TypedColumn<T> if matched, or nullptr.
   */
  template <typename T>
  TypedColumn<T> *get_column(const std::string &name) const {
    for (const auto &col : data) {
      if (col->get_name() == name && col->get_type() == typeid(T)) {
        return dynamic_cast<TypedColumn<T> *>(col.get());
      }
    }
    return nullptr;
  }

  ReflectionTable select(const std::vector<size_t> &selected_rows) const {
    ReflectionTable filtered;
    filtered.h5_filepath = this->h5_filepath;

    for (const auto &col : data) {
      filtered.data.push_back(col->clone_filtered(selected_rows));
    }

    return filtered;
  }

  template <typename T>
  ReflectionTable
  select(const std::string &column_name,
         std::function<bool(const std::span<const T>)> predicate) const {
    auto *col = get_column<T>(column_name);
    if (!col)
      throw std::runtime_error("Column not found or wrong type: " +
                               column_name);

    const auto &span = col->span();
    const auto &shape = col->get_shape();
    size_t rows = shape[0];
    size_t cols = shape.size() > 1 ? shape[1] : 1;

    std::vector<size_t> selected_rows;

    for (size_t i = 0; i < rows; ++i) {
      std::span<const T> row(&span(i, 0), cols);
      if (predicate(row)) {
        selected_rows.push_back(i);
      }
    }

    return select(selected_rows); // Reuse index-based select logic
  }

  ReflectionTable select(const std::vector<bool> &mask) const {
    std::vector<size_t> selected_rows;
    for (size_t i = 0; i < mask.size(); ++i) {
      if (mask[i]) {
        selected_rows.push_back(i);
      }
    }
    return select(selected_rows);
  }

  template <typename T>
  void add_column(const std::string &name, const size_t rows, const size_t cols,
                  const std::vector<T> &column_data) {
    add_column(name, std::vector<size_t>{rows, cols}, column_data);
  }

  template <typename T>
  void add_column(const std::string &name, const std::vector<size_t> &shape,
                  const std::vector<T> &column_data) {
    auto col = std::make_unique<TypedColumn<T>>(name, shape, column_data);

    // Ensure row count consistency
    if (!data.empty() && col->get_shape()[0] != get_row_count()) {
      throw std::runtime_error("Row count mismatch when adding column: " +
                               name);
    }

    data.push_back(std::move(col));
  }

  // void write(const std::string &filename,
  //            const std::string &group = "/dials/processing/group_0") const {
  //   for (const auto &col : data) {
  //     if (col->get_type() == typeid(double)) {
  //       auto *typed = dynamic_cast<const TypedColumn<double> *>(col.get());
  //       write_data_to_h5_file(filename, group + "/" + col->get_name(),
  //                             typed->data);
  //     } else if (col->get_type() == typeid(int)) {
  //       auto *typed = dynamic_cast<const TypedColumn<int> *>(col.get());
  //       write_data_to_h5_file(filename, group + "/" + col->get_name(),
  //                             typed->data);
  //     } else {
  //       std::cerr << "Skipping write for unsupported type: " <<
  //       col->get_name()
  //                 << "\n";
  //     }
  //   }
  // }
};
#pragma endregion