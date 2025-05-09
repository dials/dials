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
#include "dx2/h5/h5utils.hpp"
#include "dx2/h5/h5write.hpp"
#include "dx2/logging.hpp"
#include <experimental/mdspan>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace h5dispatch;
using namespace h5utils;

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

#pragma region Predicate Helpers
/**
 * @brief Represents a type-safe row predicate bound to a specific
 * column.
 *
 * This struct is used to associate a named column with a predicate
 * function that inspects values from the corresponding column. It
 * enables statically typed, variadic filtering of table rows.
 *
 * @tparam T The expected type of the column.
 */
template <typename T> struct ColumnPredicate {
  ///< The name of the column to apply the predicate on.
  std::string column_name;
  ///< A predicate function over the column mdspan and row index.
  std::function<bool(const mdspan_type<T> &, size_t)> predicate;
};

namespace logic {

/**
 * @brief Logical combination strategy for find_rows.
 */
enum class LogicalOp {
  Or, ///< Union of all matching rows (default)
  And ///< Intersection of all matching rows
};

} // namespace logic
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

  /**
   * @brief Merges a vector of row indices into a set.
   *
   * This is used to accumulate matches from multiple predicates
   * during a logical OR combination.
   *
   * @param set The target set to insert into.
   * @param rows The row indices to insert.
   */
  void merge_into_set(std::unordered_set<size_t> &set,
                      const std::vector<size_t> &rows) const {
    set.insert(rows.begin(), rows.end());
  }

  /**
   * @brief Evaluate a single ColumnPredicate against the table.
   *
   * This method retrieves the column referenced by the predicate,
   * and evaluates the predicate function for each row. A list of
   * row indices for which the predicate returned true is returned.
   *
   * The column is accessed via `column<T>(...)` which returns an
   * `mdspan` view into the column data. The predicate receives this
   * view and an index `i`, and should return true if that row is to be kept.
   *
   * Example predicate:
   *   [](const auto &span, size_t i) { return span(i, 2) > 1.0; }
   * will filter rows where the third component (z) exceeds 1.0.
   *
   * @tparam T The type of the column.
   * @param condition The column predicate to evaluate.
   * @return A vector of matching row indices.
   * @throws std::runtime_error if the column is missing or of the wrong type.
   */
  template <typename T>
  std::vector<size_t>
  evaluate_column_predicate(const ColumnPredicate<T> &condition) const {
    auto opt = column<T>(condition.column_name);
    if (!opt)
      throw std::runtime_error("Column not found or wrong type: " +
                               condition.column_name);

    const mdspan_type<T> &span = *opt;
    std::vector<size_t> matches;
    for (size_t i = 0; i < span.extent(0); ++i) {
      if (condition.predicate(span, i)) {
        matches.push_back(i);
      }
    }
    return matches;
  }

  /**
   * @brief Type trait used to ensure all variadic arguments are
   * ColumnPredicate<T>.
   */
  template <typename T> struct is_column_predicate : std::false_type {};

  template <typename T>
  struct is_column_predicate<ColumnPredicate<T>> : std::true_type {};

public:
#pragma region Constructors
  ReflectionTable() = default;

  /**
   * @brief Constructs a ReflectionTable with given experiment IDs and
   * identifiers.
   *
   * @param experiment_ids A vector of experiment IDs.
   * @param identifiers A vector of identifiers.
   */
  ReflectionTable(const std::vector<uint64_t> &experiment_ids,
                  const std::vector<std::string> &identifiers)
      : experiment_ids(experiment_ids), identifiers(identifiers) {}

  /**
   * @brief Constructs a ReflectionTable from an HDF5 file.
   *
   * This constructor loads datasets from the specified HDF5 file and
   * populates the table with the data.
   *
   * @param h5_filepath The path to the HDF5 file.
   */
  ReflectionTable(const std::string &h5_filepath) : h5_filepath(h5_filepath) {
    auto start = std::chrono::high_resolution_clock::now(); // ⏱ Start timer

    // Discover all datasets in the default reflection group
    std::vector<std::string> datasets =
        get_datasets_in_group(h5_filepath, DEFAULT_REFL_GROUP);

    if (datasets.empty()) {
      dx2_log::warning(
          fmt::format("No datasets found in group '{}'", DEFAULT_REFL_GROUP));
    }

    // Open the HDF5 file
    H5File file(H5Fopen(h5_filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (!file) {
      throw std::runtime_error("Could not open file: " + h5_filepath);
    }

    // Open group and read experiment metadata
    H5Group group(H5Gopen2(file, DEFAULT_REFL_GROUP.c_str(), H5P_DEFAULT));
    if (group) {
      read_experiment_metadata(group, experiment_ids, identifiers);
    }

    // Loop over every dataset path in the group
    for (const auto &dataset : datasets) {
      std::string dataset_name = get_dataset_name(dataset);

      // Open the specific dataset (within the file opened above)
      H5Dataset dataset_id(H5Dopen2(file, dataset.c_str(), H5P_DEFAULT));
      if (!dataset_id) {
        dx2_log::warning(fmt::format("Could not open dataset '{}'", dataset));
        continue;
      }

      try {
        h5dispatch::dispatch_h5_dataset_type(dataset_id, [&](auto tag) {
          using T = typename decltype(tag)::type;
          auto result =
              read_array_with_shape_from_h5_file<T>(h5_filepath, dataset);

          data.push_back(std::make_unique<TypedColumn<T>>(
              dataset_name, result.shape, result.data));

          dx2_log::debug("Loaded column: {} with type {}", dataset_name,
                         typeid(T).name());
        });
      } catch (const std::exception &e) {
        dx2_log::warning(
            fmt::format("Skipping dataset '{}': {}", dataset, e.what()));
        // Continue to the next dataset
      }
    }

    dx2_log::debug("Loaded {} column(s) from group '{}'", data.size(),
                   DEFAULT_REFL_GROUP);

    auto end = std::chrono::high_resolution_clock::now(); // ⏱ End timer
    dx2_log::debug("ReflectionTable loaded in {:.4f}s",
                   std::chrono::duration<double>(end - start).count());
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
   * @brief Set the list of experiment IDs.
   */
  void set_experiment_ids(const std::vector<uint64_t> &ids) {
    experiment_ids = ids;
  }

  /**
   * @brief Get the list of identifiers.
   */
  const std::vector<std::string> &get_identifiers() const {
    return identifiers;
  }

  /**
   * @brief Set the list of identifiers.
   */
  void set_identifiers(const std::vector<std::string> &ids) {
    identifiers = ids;
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
  std::optional<mdspan_type<T>> column(const std::string &name) const {
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
   * @brief Finds all rows matching any or all of the given typed
   * predicates.
   *
   * This function accepts any number of typed predicates, each bound to
   * a specific column. It returns the union or intersection of all
   * matching rows, depending on the specified logical operation.
   *
   * Uses a parameter pack to accept a variadic number of predicates,
   * and expands them into individual calls to evaluate_column_predicate.
   * The resulting vectors are combined based on the logical operator.
   *
   * @tparam Predicates Parameter pack of ColumnPredicate<T> types.
   * @param op Logical combination (AND or OR).
   * @param conds The predicates to apply to rows.
   * @return A vector of row indices matching the logical combination.
   */
  template <typename... Predicates>
  std::vector<size_t> find_rows(logic::LogicalOp op,
                                const Predicates &...conds) const {
    // Ensure all predicates are of the correct type
    static_assert((is_column_predicate<Predicates>::value && ...),
                  "All arguments to find_rows must be ColumnPredicate<T>");

    // Expand the parameter pack into a vector of row-index vectors
    std::vector<std::vector<size_t>> all_matches = {
        evaluate_column_predicate(conds)...};
    std::unordered_set<size_t> result_set;

    if (op == logic::LogicalOp::Or) {
      // Combine all matching row indices (set union)
      for (const auto &matches : all_matches) {
        result_set.insert(matches.begin(), matches.end());
      }
    } else if (op == logic::LogicalOp::And) {
      // Combine only common row indices (set intersection)
      if (all_matches.empty())
        return {};

      // Start with the first match set as the initial intersection base
      std::unordered_set<size_t> intersection(all_matches[0].begin(),
                                              all_matches[0].end());

      // Iteratively refine the intersection by retaining only elements
      // that are also present in the next predicate's results
      for (size_t i = 1; i < all_matches.size(); ++i) {
        std::unordered_set<size_t> next(all_matches[i].begin(),
                                        all_matches[i].end());

        // Iterate over current intersection set
        for (auto it = intersection.begin(); it != intersection.end();) {
          // If the current value is not present in the next result set,
          // erase it from the intersection set
          if (next.find(*it) == next.end()) {
            it = intersection.erase(it); // erase returns the next iterator
          } else {
            ++it; // move to the next element if it's still present
          }
        }
      }

      result_set = std::move(intersection);
    }

    return std::vector<size_t>(result_set.begin(), result_set.end());
  }

  /**
   * @brief Finds all rows matching all of the given typed predicates.
   *
   * This overload enables simpler syntax when AND is desired (default).
   * Because default parameters don't work with parameter packs,
   * we provide this wrapper to make the interface intuitive.
   *
   * @tparam Predicates Parameter pack of ColumnPredicate<T> types.
   * @param conds The predicates to apply to rows.
   * @return A vector of row indices matching all predicates.
   */
  template <typename... Predicates>
  std::vector<size_t> find_rows(const Predicates &...conds) const {
    return find_rows(logic::LogicalOp::And, conds...);
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
   * @brief Adds a new column to the table with arbitrary shape.
   *
   * This is the most general column-adding interface. You must specify
   * the full shape of the column explicitly. The column is stored using
   * a flat data buffer matching this shape in row-major order.
   *
   * @tparam T The data type of the column (e.g., `double`, `int`,
   * `std::array<T, N>`).
   * @param name The name of the column.
   * @param shape A vector describing the shape of the column (e.g.,
   * `{N}` for 1D, `{N, M}` for 2D).
   * @param column_data A flat vector of `T` values, whose total number
   * of elements must match the product of `shape`.
   *
   * @throws std::runtime_error if shape does not match existing row
   * count (if table is non-empty).
   */
  template <typename T>
  void add_column(const std::string &name, const std::vector<size_t> &shape,
                  const std::vector<T> &column_data) {
    auto col = std::make_unique<TypedColumn<T>>(name, shape, column_data);

    // Check for duplicate column names
    for (const auto &existing_col : data) {
      if (existing_col->get_name() == name) {
        throw std::runtime_error("Column with name already exists: " + name);
      }
    }

    // Check if type T is supported
    const auto &registry = h5dispatch::get_supported_types();
    bool supported = std::any_of(registry.begin(), registry.end(),
                                 [](const h5dispatch::H5TypeInfo &info) {
                                   return info.cpp_type == typeid(T);
                                 });

    if (!supported) {
      throw std::runtime_error(
          "Attempted to add column with unsupported type: " +
          std::string(typeid(T).name()));
    }

    // Ensure row count consistency
    if (!data.empty() && col->get_shape()[0] != get_row_count()) {
      throw std::runtime_error("Row count mismatch when adding column: " +
                               name);
    }

    // Add the new column to the table
    data.push_back(std::move(col));
  }

  /**
   * @brief Adds a 1D column to the table.
   *
   * This is a convenience overload for scalar or vector-style columns
   * where each row has a single value. Internally calls the shape-based
   * overload with a shape of `{N}`.
   *
   * @tparam T The data type of the column (e.g., `int`, `double`).
   * @param name The name of the column.
   * @param column_data A vector of `T` values, one per row.
   */
  template <typename T>
  void add_column(const std::string &name, const std::vector<T> &column_data) {
    add_column(name, std::vector<size_t>{column_data.size()}, column_data);
  }

  /**
   * @brief Adds a 2D column to the table from row and column dimensions.
   *
   * Use this for matrix-style data where each row is a fixed-length vector.
   * Internally calls the shape-based overload with a shape of `{rows, cols}`.
   *
   * @tparam T The data type of the column elements.
   * @param name The name of the column.
   * @param rows Number of rows (typically must match existing row count).
   * @param cols Number of columns per row.
   * @param column_data Flat vector of values in row-major order (size must be
   * `rows * cols`).
   */
  template <typename T>
  void add_column(const std::string &name, const size_t rows, const size_t cols,
                  const std::vector<T> &column_data) {
    add_column(name, std::vector<size_t>{rows, cols}, column_data);
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
  void write(std::string_view filename,
             std::string_view group = "/dials/processing/group_0") const {
    std::string fname(filename);
    std::string gpath(group);

    // Suppress errors when opening non-existent files, groups, datasets..
    H5ErrorSilencer silencer;

    // 🗂️ Ensure the file exists or create it before writing
    H5File file(H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    if (!file) {
      file = H5File(
          H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      if (!file) {
        throw std::runtime_error("Failed to create or open file: " + fname);
      }
    }

    // Open or create group
    H5Group group_id = traverse_or_create_groups(file, gpath);
    if (!group_id) {
      throw std::runtime_error("Failed to create or open group: " + gpath);
    }

    // Write metadata
    write_experiment_metadata(group_id, experiment_ids, identifiers);

    // 🔁 Write all columns
    for (const auto &col : data) {
      // 🏗️ Construct full dataset path: group + column name
      const std::string &name = col->get_name();

      // Define a lambda that writes a column of a specific type T
      auto write_col = [&](auto tag) {
        using T = typename decltype(tag)::type;

        // 🧪 Try to cast the typeless base pointer to TypedColumn<T>
        const auto *typed = col->as<T>();
        if (!typed) {
          // This should not happen unless type registry is inconsistent
          dx2_log::error(
              fmt::format("Internal type mismatch for column '{}'", name));
          return;
        }

        // 💾 Call the HDF5 writer to write raw data to file
        write_raw_data_to_h5_group<T>(
            group_id, name, typed->data.data(),
            std::vector<hsize_t>(typed->shape.begin(), typed->shape.end()));
      };

      // 🌀 Dispatch the column type and invoke the write_col lambda
      try {
        dispatch_column_type(col->get_type(), write_col);
      } catch (const std::exception &e) {
        // ⚠️ If the type is unsupported or an error occurs, print warning
        dx2_log::warning(fmt::format("Skipping column {}: {}", name, e.what()));
      }
    }
  }
#pragma endregion
};
#pragma endregion