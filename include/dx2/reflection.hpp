/**
 * @file reflection.hpp
 * @brief Lightweight, type-safe reflection table for HDF5 datasets.
 *
 * This header defines a flexible structure for interacting with HDF5
 * datasets as typed reflection tables. It allows datasets to be loaded
 * into typed polymorphic containers, filtered by row selection, queried
 * by column, and written back to disk.
 *
 * - `ReflectionTable` manages a collection of columns loaded from an
 *   HDF5 group.
 * - Each column is stored as a type-erased `ColumnBase`, which can be
 *   downcast to a specific `TypedColumn<T>` when needed.
 * - Support for adding, querying, selecting, and writing columns.
 */
#pragma once

#include "dx2/h5/h5dispatch.hpp"
#include "dx2/h5/h5read_processed.hpp"
#include "dx2/h5/h5write.hpp"
#include <experimental/mdspan>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace h5dispatch;

#pragma region Type Helpers
/*
 * Type aliases for centralised control of storage types.
 */

// Internal storage type for the data
template <typename T> using storage_type = std::vector<T>;
// Type alias for mdspan (currently using backport std::experimental::mdspan)
template <typename T>
using mdspan_type =
    std::experimental::mdspan<T, std::experimental::dextents<size_t, 2>>;
#pragma endregion

#pragma region Column Definition
// Forward declaration of TypedColumn
template <typename T> struct TypedColumn;

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

  /**
   * @brief Attempts to cast this column to a TypedColumn<T>.
   *
   * @tparam T The target type to cast to.
   * @return A pointer to the TypedColumn<T>, or nullptr if types don't match.
   */
  template <typename T> const TypedColumn<T> *as() const {
    if (get_type() != std::type_index(typeid(T)))
      return nullptr;
    return dynamic_cast<const TypedColumn<T> *>(this);
  }
};

/**
 * @brief Specialized container for HDF5 data of a specific type.
 */
template <typename T> struct TypedColumn : public ColumnBase {
  /**
   * @brief Constructs a column with a given name and shape, allocating storage.
   *
   * @param name The name of the column.
   * @param shape A vector specifying the shape (e.g. [rows, cols]).
   */
  TypedColumn(const std::string &name, const std::vector<size_t> &shape)
      : name(name), shape(shape), data(compute_size(shape)) {
    setup_mdspan();
  }

  /**
   * @brief Constructs a column with a given name, shape, and preloaded data.
   *
   * @param name The name of the column.
   * @param shape Shape of the column.
   * @param data Flat vector of values, matching the given shape.
   */
  TypedColumn(const std::string &name, const std::vector<size_t> &shape,
              const std::vector<T> &data)
      : name(name), shape(shape), data(data) {
    setup_mdspan();
  }

  /**
   * @brief Constructs a column from row/col dimensions.
   *
   * @param name The column name.
   * @param rows Number of rows.
   * @param cols Number of columns.
   */
  TypedColumn(const std::string &name, size_t rows, size_t cols)
      : name(name), shape{rows, cols}, data(rows * cols) {
    setup_mdspan();
  }

  std::vector<size_t> get_shape() const override { return shape; }
  std::string get_name() const override { return name; }
  std::type_index get_type() const override {
    return std::type_index(typeid(T));
  }

  /**
   * @brief Returns an mdspan view over the internal data.
   *
   * @return A constant reference to the mdspan view of the data.
   */
  const mdspan_type<T> &span() const { return mdspan_data; }

  /**
   * @brief Creates a filtered copy of the column, keeping only selected rows.
   *
   * This method returns a new TypedColumn containing only the rows specified
   * by the input indices. Column structure and type are preserved.
   *
   * @param selected_rows Indices of rows to keep.
   * @return A new column with only the specified rows.
   */
  std::unique_ptr<ColumnBase>
  clone_filtered(const std::vector<size_t> &selected_rows) const override {
    size_t new_rows = selected_rows.size();
    size_t cols = shape.size() > 1 ? shape[1] : 1;

    std::vector<T> new_data;
    new_data.reserve(new_rows * cols);

    for (size_t row : selected_rows) {
      for (size_t col = 0; col < cols; ++col) {
        // 🧵 Extract value row-wise
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
  /**
   * @brief Initializes the mdspan view based on the current shape.
   */
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

  /**
   * @brief Computes the total number of elements for a given shape.
   *
   * @param shape A vector of dimensions.
   * @return The product of all shape dimensions.
   */
  static size_t compute_size(const std::vector<size_t> &shape) {
    if (shape.empty())
      return 0;
    size_t s = 1;
    for (auto dim : shape)
      s *= dim;
    return s;
  }
};

template <typename Callback>
void dispatch_column_type(const ColumnBase &col, Callback &&cb) {
  dispatch_column_type(col.get_type(), std::forward<Callback>(cb));
}
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
  std::vector<uint64_t> experiment_ids;
  std::vector<std::string> identifiers;

  /// Get the number of rows in the first column, assuming all columns
  /// have the same number of rows
  size_t get_row_count() const {
    if (data.empty())
      return 0;
    return data.front()->get_shape()[0];
  }

  /// Check if a column with the given name and type exists
  template <typename T>
  bool matches_column(const ColumnBase &col, const std::string &name) const {
    return col.get_name() == name && col.get_type() == typeid(T);
  }

public:
#pragma region Constructors
  ReflectionTable() = default;

  ReflectionTable(const std::string &h5_filepath) : h5_filepath(h5_filepath) {
    // Discover all datasets in the default reflection group
    std::vector<std::string> datasets =
        get_datasets_in_group(h5_filepath, DEFAULT_REFL_GROUP);

    // Open the HDF5 file
    hid_t file = H5Fopen(h5_filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      throw std::runtime_error("Could not open file: " + h5_filepath);
    }

    // Open group and read experiment metadata
    hid_t group = H5Gopen(file, DEFAULT_REFL_GROUP.c_str(), H5P_DEFAULT);
    if (group >= 0) {
      read_experiment_metadata(group, experiment_ids, identifiers);
      H5Gclose(group);
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
        h5dispatch::dispatch_h5_dataset_type(dataset_id, [&](auto tag) {
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
#pragma endregion

#pragma region Metadata Access
  /**
   * @brief Get the list of experiment IDs.
   */
  const std::vector<uint64_t> &get_experiment_ids() const {
    return experiment_ids;
  }

  /**
   * @brief Get the list of identifiers.
   */
  const std::vector<std::string> &get_identifiers() const {
    return identifiers;
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
#pragma endregion

#pragma region Column Access
  /**
   * @brief Retrieves a read-only view (mdspan) of a typed column by
   * name.
   *
   * This function attempts to locate a column in the table that matches
   * both the given name and expected type `T`. If a matching column is
   * found, it returns a non-owning reference wrapper around the
   * `mdspan` view of that column’s data. If no matching column is
   * found, it returns std::nullopt.
   *
   * @tparam T The expected type of the column (e.g., int, double).
   * @param name The name of the column to retrieve.
   * @return An optional containing a reference to the mdspan view, or
   * std::nullopt if not found.
   *
   * @note The returned mdspan reference is only valid as long as the
   * owning ReflectionTable instance and its internal data remain alive.
   */
  template <typename T>
  std::optional<mdspan_type<T>> get_column(const std::string &name) const {
    // Iterate though all columns in the table
    for (const auto &col : data) {
      // Check if the column matches the requested type and name
      if (matches_column<T>(*col, name)) {
        // Attempt to cast to the specific type
        auto *typed = dynamic_cast<TypedColumn<T> *>(col.get());
        if (typed) { // If cast is successful, return the span
          return typed->span();
        }
      }
    }
    // If no matching column is found, return an empty optional
    return std::nullopt;
  }

  /**
   * @brief Finds rows matching a predicate.
   *
   * Applies the given predicate to each row (represented by its index).
   * The predicate can fetch any columns it needs via the table interface.
   *
   * @param predicate A function taking a row index and returning a bool.
   * @return A vector of row indices for which the predicate returned true.
   */
  std::vector<size_t> find_rows(std::function<bool(size_t)> predicate) const {
    std::vector<size_t> matching_rows;
    size_t row_count = get_row_count();
    for (size_t i = 0; i < row_count; ++i) {
      if (predicate(i)) {
        matching_rows.push_back(i);
      }
    }
    return matching_rows;
  }

  /**
   * @brief Returns a new ReflectionTable with only the specified rows.
   *
   * @param selected_rows The list of row indices to retain.
   * @return A new ReflectionTable containing only the selected rows.
   */
  ReflectionTable select(const std::vector<size_t> &selected_rows) const {
    ReflectionTable filtered;
    filtered.h5_filepath = this->h5_filepath;

    for (const auto &col : data) {
      filtered.data.push_back(col->clone_filtered(selected_rows));
    }

    // Copy experiment_ids and identifiers
    filtered.experiment_ids = this->experiment_ids;
    filtered.identifiers = this->identifiers;

    return filtered;
  }

  /**
   * @brief Returns a new ReflectionTable from a boolean mask.
   *
   * Converts a boolean mask into row indices and calls index-based select.
   *
   * @param mask A boolean mask, where `true` indicates a row to keep.
   * @return A new ReflectionTable containing only the selected rows.
   */
  ReflectionTable select(const std::vector<bool> &mask) const {
    std::vector<size_t> selected_rows;
    for (size_t i = 0; i < mask.size(); ++i) {
      if (mask[i]) {
        selected_rows.push_back(i);
      }
    }
    return select(selected_rows);
  }
#pragma endregion

#pragma region Column Modification
  /**
   * @brief Adds a new column to the table from row/col dimensions and flat
   * data.
   *
   * @tparam T The data type of the column.
   * @param name The name of the column.
   * @param rows Number of rows.
   * @param cols Number of columns.
   * @param column_data Flat vector of data values.
   */
  template <typename T>
  void add_column(const std::string &name, const size_t rows, const size_t cols,
                  const std::vector<T> &column_data) {
    add_column(name, std::vector<size_t>{rows, cols}, column_data);
  }

  /**
   * @brief Adds a new column to the table from shape vector and flat data.
   *
   * @tparam T The data type of the column.
   * @param name The name of the column.
   * @param shape Shape of the column.
   * @param column_data Flat vector of data values.
   *
   * @throws std::runtime_error if shape does not match existing row count.
   */
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
#pragma endregion

#pragma region Write
  /**
   * @brief Writes all columns to an HDF5 file under the given group.
   *
   * This function iterates through all stored columns and writes each one
   * to disk, provided its type is supported by the HDF5 backend.
   * Currently supported types include: double, int, int64_t, and uint64_t.
   *
   * Columns with unsupported types will be skipped with a warning.
   *
   * @param filename Output HDF5 file path.
   * @param group Target group inside the HDF5 file (defaults to
   * /dials/processing/group_0).
   */
  void write(const std::string &filename,
             const std::string &group = "/dials/processing/group_0") const {
    // 🗂️ Ensure the file exists or create it before writing
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
      file =
          H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (file < 0) {
        throw std::runtime_error("Failed to create or open file: " + filename);
      }
    }

    // Traverse and get group
    hid_t group_id = traverse_or_create_groups(file, group);

    // Write metadata
    write_experiment_metadata(group_id, experiment_ids, identifiers);

    H5Gclose(group_id);

    // 🔁 Write all data columns
    for (const auto &col : data) {
      // 🏗️ Construct full dataset path: group + column name
      std::string path = group + "/" + col->get_name();

      // Define a lambda that writes a column of a specific type T
      auto write_col = [&](auto tag) {
        using T = typename decltype(tag)::type;

        // 🧪 Try to cast the typeless base pointer to TypedColumn<T>
        const auto *typed = col->as<T>();
        if (!typed) {
          // This should not happen unless type registry is inconsistent
          std::cerr << "Internal type mismatch for column: " << col->get_name()
                    << "\n";
          return;
        }

        // 💾 Call the HDF5 writer to write raw data to file
        write_raw_data_to_h5_file<T>(
            filename, path, typed->data.data(),
            {typed->shape.begin(), typed->shape.end()});
      };

      // 🌀 Dispatch the column type and invoke the write_col lambda
      try {
        dispatch_column_type(col->get_type(), write_col);
      } catch (const std::exception &e) {
        // ⚠️ If the type is unsupported or an error occurs, print warning
        std::cerr << "Skipping column " << col->get_name() << ": " << e.what()
                  << "\n";
      }
    }

    H5Fclose(file);
  }
#pragma endregion
};
#pragma endregion