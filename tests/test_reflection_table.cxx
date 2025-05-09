#include <dx2/reflection.hpp>
#include <filesystem>
#include <gtest/gtest.h>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

/*
 * Test fixture for ReflectionTable tests. This fixture sets up the test
 * environment for the ReflectionTable class and provides a common path
 * to the test file.
 */
class ReflectionTableTest : public ::testing::Test {
protected:
  std::filesystem::path test_file_path;

  void SetUp() override {
    test_file_path = std::filesystem::current_path() / "data/cut_strong.refl";
  }
};

#pragma region Loading
TEST_F(ReflectionTableTest, LoadsAllColumns) {
  ReflectionTable table(test_file_path.string());

  auto names = table.get_column_names();
  std::cout << "Loaded column names:\n";
  for (const auto &name : names) {
    std::cout << "  - " << name << "\n";
  }

  EXPECT_FALSE(names.empty());
  EXPECT_NE(std::find(names.begin(), names.end(), "xyzobs.px.value"),
            names.end());
}

TEST_F(ReflectionTableTest, TryGetColumnSucceedsOnCorrectType) {
  ReflectionTable table(test_file_path.string());

  auto px = table.column<double>("xyzobs.px.value");
  ASSERT_TRUE(px.has_value());

  const auto &span = px.value();
  std::cout << "xyzobs.px.value shape: " << span.extent(0) << "x"
            << span.extent(1) << "\n";
  EXPECT_EQ(span.extent(1), 3);
  EXPECT_GT(span.extent(0), 0);
}

TEST_F(ReflectionTableTest, TryGetColumnFailsOnWrongType) {
  ReflectionTable table(test_file_path.string());

  auto col = table.column<int>("xyzobs.px.value");
  std::cout << "Column retrieval (int): "
            << (col.has_value() ? "found" : "not found") << "\n";
  EXPECT_FALSE(col.has_value());
}

TEST_F(ReflectionTableTest, ConstructorWithIDsStoresMetadataCorrectly) {
  std::vector<uint64_t> ids = {123, 456};
  std::vector<std::string> labels = {"foo", "bar"};

  ReflectionTable table(ids, labels);

  EXPECT_EQ(table.get_experiment_ids(), ids);
  EXPECT_EQ(table.get_identifiers(), labels);
}
#pragma endregion

#pragma region Validation
TEST_F(ReflectionTableTest, AllDatasetsLoadSuccessfully) {
  ReflectionTable table(test_file_path.string());

  std::vector<std::string> expected = get_datasets_in_group(
      test_file_path.string(), "/dials/processing/group_0");

  std::vector<std::string> loaded = table.get_column_names();

  for (const std::string &full_dataset_path : expected) {
    std::string name = get_dataset_name(full_dataset_path);

    // Make sure it's present
    auto it = std::find(loaded.begin(), loaded.end(), name);
    EXPECT_NE(it, loaded.end()) << "Dataset not loaded: " << name;

    if (it != loaded.end()) {
      std::cout << "✅ Dataset loaded: " << name << "\n";
    } else {
      std::cerr << "❌ Dataset missing: " << name << "\n";
    }
  }
}
#pragma endregion

#pragma region Access
TEST_F(ReflectionTableTest, AccessNonExistentColumn) {
  ReflectionTable table(test_file_path.string());

  // Attempt to access a column that does not exist
  auto non_existent_column = table.column<double>("non_existent_column");

  /**
   * Attemtping to access a non-existent column should return an empty
   * optional.
   *
   * This is a test to ensure that the column retrieval mechanism
   * correctly handles cases where the requested column does not exist.
   */

  std::cout << "Attempting to access non-existent column: "
            << (non_existent_column.has_value() ? "found" : "not found")
            << "\n";

  // Ensure the column retrieval fails
  EXPECT_FALSE(non_existent_column.has_value());
}
#pragma endregion

#pragma region Adding
TEST_F(ReflectionTableTest, AddColumn1D) {
  ReflectionTable table;
  std::vector<double> data{1.0, 2.0, 3.0};
  table.add_column("col_1d", data); // 1D: implicit shape [3]

  auto col = table.column<double>("col_1d");
  ASSERT_TRUE(col.has_value());

  std::cout << "[AddColumn1D] Shape: " << col->extent(0) << "x"
            << col->extent(1) << "\n";
  for (size_t i = 0; i < col->extent(0); ++i) {
    std::cout << "  col(" << i << ", 0) = " << (*col)(i, 0) << "\n";
  }
}

TEST_F(ReflectionTableTest, AddColumn2D) {
  ReflectionTable table;
  std::vector<double> data{1.1, 1.2, 2.1, 2.2, 3.1, 3.2}; // 2D: shape [3, 2]
  table.add_column("col_2d", 3, 2, data);

  auto col = table.column<double>("col_2d");
  ASSERT_TRUE(col.has_value());

  std::cout << "[AddColumn2D] Shape: " << col->extent(0) << "x"
            << col->extent(1) << "\n";
  for (size_t i = 0; i < col->extent(0); ++i) {
    for (size_t j = 0; j < col->extent(1); ++j) {
      std::cout << "  col(" << i << ", " << j << ") = " << (*col)(i, j) << "\n";
    }
  }
}

TEST_F(ReflectionTableTest, AddColumnExplicit) {
  ReflectionTable table;

  // 2D shape explicitly specified as {3, 3}
  std::vector<double> data{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  table.add_column("col_explicit", std::vector<size_t>{3, 3}, data);

  auto col = table.column<double>("col_explicit");
  ASSERT_TRUE(col.has_value());

  std::cout << "[AddColumnExplicit] Shape: " << col->extent(0) << "x"
            << col->extent(1) << "\n";
  for (size_t i = 0; i < col->extent(0); ++i) {
    for (size_t j = 0; j < col->extent(1); ++j) {
      std::cout << "  col(" << i << ", " << j << ") = " << (*col)(i, j) << "\n";
    }
  }

  // Basic value check
  EXPECT_DOUBLE_EQ((*col)(1, 1), 5.0);
  EXPECT_DOUBLE_EQ((*col)(2, 2), 9.0);
}

TEST_F(ReflectionTableTest, AddColumnRejectsRank3) {
  ReflectionTable table;

  // Invalid rank-3 shape: {2, 2, 3}
  std::vector<double> data(12);
  std::iota(data.begin(), data.end(), 0.0);

  try {
    table.add_column("col_rank3", std::vector<size_t>{2, 2, 3}, data);
    FAIL() << "Expected std::runtime_error";
  } catch (const std::runtime_error &e) {
    std::cout << "[AddColumnRejectsRank3] Caught error: " << e.what() << "\n";
    EXPECT_TRUE(std::string(e.what()).find("unsupported rank") !=
                std::string::npos);
  } catch (...) {
    FAIL() << "Expected std::runtime_error: unsupported rank";
  }
}

TEST_F(ReflectionTableTest, AddColumnRejectsMismatchedRowCount) {
  ReflectionTable table;

  // First column: 3 rows
  table.add_column("col1", std::vector<double>{1.0, 2.0, 3.0});

  // Second column: 2 rows – should throw
  std::vector<double> bad_data{4.0, 5.0};
  try {
    table.add_column("col2", bad_data);
    FAIL() << "Expected std::runtime_error due to row count mismatch";
  } catch (const std::runtime_error &e) {
    std::cout << "[AddColumnRejectsMismatchedRowCount] Caught error: "
              << e.what() << "\n";
    EXPECT_TRUE(std::string(e.what()).find("Row count mismatch") !=
                std::string::npos);
  } catch (...) {
    FAIL() << "Expected std::runtime_error: row count mismatch";
  }
}

TEST_F(ReflectionTableTest, AddDuplicateColumnThrows) {
  ReflectionTable table;
  std::vector<double> data{1.0, 2.0, 3.0};
  table.add_column("dup", data);

  try {
    table.add_column("dup", data); // Add column with same name again
    FAIL() << "Expected logic_error or appropriate overwrite warning";
  } catch (const std::exception &e) {
    std::cout << "[AddDuplicateColumnThrows] Caught: " << e.what() << "\n";
    EXPECT_TRUE(std::string(e.what()).find("already exists") !=
                std::string::npos);
  }
}

TEST_F(ReflectionTableTest, AddUnsupportedColumnTypeThrows) {
  ReflectionTable table;

  // Define a dummy type that is not supported
  struct DummyType {
    int x;
  };

  // Attempt to add a column with the unsupported type
  std::vector<DummyType> data(3, DummyType{42});
  try {
    table.add_column("dummy", std::vector<size_t>{3}, data);
  } catch (const std::exception &e) {
    std::cout << "[AddUnsupportedColumnTypeThrows] Error (expected): "
              << e.what() << "\n";
    SUCCEED();
    return;
  }
  FAIL() << "Expected std::runtime_error: unsupported column type";
}

TEST_F(ReflectionTableTest, WriteTableFromScratchAndReload) {
  // Construct a table from scratch
  ReflectionTable table;

  // Add an integer column
  std::vector<int> example_data{0, 2, 3, 1};
  std::string column_name = "id";
  table.add_column<int>(column_name, 4, 1, example_data);

  // Add metadata
  std::vector<uint64_t> experiment_ids{101, 202};
  std::vector<std::string> identifiers{"alpha", "beta"};
  table.set_experiment_ids(experiment_ids);
  table.set_identifiers(identifiers);

  // Write to a temporary file
  std::filesystem::path temp_file =
      std::filesystem::current_path() / "reflection_test_scratch_write.h5";
  table.write(temp_file.string());

  // Reload from file
  ReflectionTable loaded(temp_file.string());

  // Validate column presence and contents
  auto span = loaded.column<int>("id");
  ASSERT_TRUE(span.has_value());
  ASSERT_EQ(span->extent(0), 4);
  ASSERT_EQ(span->extent(1), 1);
  for (size_t i = 0; i < 4; ++i) {
    EXPECT_EQ((*span)(i, 0), example_data[i]);
  }

  // Validate metadata round-trip
  const auto &reloaded_ids = loaded.get_experiment_ids();
  const auto &reloaded_identifiers = loaded.get_identifiers();

  ASSERT_EQ(reloaded_ids.size(), experiment_ids.size());
  ASSERT_EQ(reloaded_identifiers.size(), identifiers.size());

  for (size_t i = 0; i < experiment_ids.size(); ++i) {
    EXPECT_EQ(reloaded_ids[i], experiment_ids[i]);
    EXPECT_EQ(reloaded_identifiers[i], identifiers[i]);
  }

  // Clean up
  std::filesystem::remove(temp_file);
}
#pragma endregion

#pragma region Selection
TEST_F(ReflectionTableTest, SelectWithMaskReturnsSameResultAsExplicitIndices) {
  ReflectionTable table(test_file_path.string());

  auto px = table.column<double>("xyzobs.px.value");
  ASSERT_TRUE(px);
  const auto &span = px.value();

  std::vector<bool> mask(span.extent(0), false);
  std::vector<size_t> indices;

  for (size_t i = 0; i < span.extent(0); ++i) {
    if (span(i, 2) > 1.0) {
      mask[i] = true;
      indices.push_back(i);
    }
  }

  std::cout << "Mask selected " << indices.size() << " rows\n";

  auto table_from_mask = table.select(mask);
  auto table_from_indices = table.select(indices);

  auto mask_span = table_from_mask.column<double>("xyzobs.px.value").value();
  auto idx_span = table_from_indices.column<double>("xyzobs.px.value").value();

  ASSERT_EQ(mask_span.extent(0), idx_span.extent(0));
  ASSERT_EQ(mask_span.extent(1), idx_span.extent(1));

  for (size_t i = 0; i < std::min<size_t>(mask_span.extent(0), 3); ++i) {
    std::cout << "Row " << i << ": "
              << "mask_z = " << mask_span(i, 2)
              << ", idx_z = " << idx_span(i, 2) << "\n";
  }

  for (size_t i = 0; i < mask_span.extent(0); ++i) {
    for (size_t j = 0; j < mask_span.extent(1); ++j) {
      EXPECT_EQ(mask_span(i, j), idx_span(i, j));
    }
  }
}

TEST_F(ReflectionTableTest, SelectEmptyReturnsEmptyTable) {
  ReflectionTable table(test_file_path.string());
  ReflectionTable empty = table.select(std::vector<size_t>{});

  std::cout << "Selecting with empty row set...\n";
  for (const auto &name : table.get_column_names()) {
    auto span = empty.column<double>(name);
    if (span) {
      std::cout << "  " << name << " has " << span->extent(0) << " rows\n";
      EXPECT_EQ(span.value().extent(0), 0);
    }
  }
}

TEST_F(ReflectionTableTest, SelectSubsetUsingTypedFindRows) {
  ReflectionTable table(test_file_path.string());

  // Define a predicate that selects rows where the z-coordinate (3rd column) is
  // greater than 1.0
  auto selected = table.find_rows(ColumnPredicate<double>{
      "xyzobs.px.value",
      [](const auto &span, size_t i) { return span(i, 2) > 1.0; }});

  std::cout << "Selected " << selected.size() << " rows where z > 1.0\n";
  EXPECT_FALSE(selected.empty());

  // Filter the table to include only the selected rows
  auto filtered = table.select(selected);

  // Extract the filtered xyzobs.px.value column
  auto filtered_px = filtered.column<double>("xyzobs.px.value");
  ASSERT_TRUE(filtered_px);
  const auto &filtered_span = filtered_px.value();

  // Print the first few z-values (i.e. the value at column index 2)
  for (size_t i = 0; i < std::min<size_t>(filtered_span.extent(0), 5); ++i) {
    std::cout << "z[" << i << "] = " << filtered_span(i, 2) << "\n";
  }

  // Validate all selected rows meet the predicate condition
  for (size_t i = 0; i < filtered_span.extent(0); ++i) {
    EXPECT_GT(filtered_span(i, 2), 1.0);
  }
}

TEST_F(ReflectionTableTest, SelectSubsetWithMultipleTypedPredicates_OR) {
  ReflectionTable table(test_file_path.string());

  // Combine two predicates: x > 200 (column 0) or z < 5.0 (column 2)
  auto selected =
      table.find_rows(logic::LogicalOp::Or,
                      ColumnPredicate<double>{"xyzobs.px.value",
                                              [](const auto &span, size_t i) {
                                                return span(i, 0) > 200.0;
                                              }},
                      ColumnPredicate<double>{"xyzobs.px.value",
                                              [](const auto &span, size_t i) {
                                                return span(i, 2) < 5.0;
                                              }});

  std::cout << "Selected " << selected.size()
            << " rows with x > 200 OR z < 5.0\n";
  EXPECT_FALSE(selected.empty());

  // Select rows from table and extract the column
  auto filtered = table.select(selected);
  auto filtered_px = filtered.column<double>("xyzobs.px.value");
  ASSERT_TRUE(filtered_px);
  const auto &span = filtered_px.value();

  // Verify each selected row satisfies at least one predicate
  for (size_t i = 0; i < span.extent(0); ++i) {
    EXPECT_TRUE(span(i, 0) > 200.0 || span(i, 2) < 5.0);
  }
}

TEST_F(ReflectionTableTest, SelectSubsetWithMultipleTypedPredicates_AND) {
  ReflectionTable table(test_file_path.string());

  // Combine two predicates: x > 200 (column 0) and z < 5.0 (column 2)
  auto selected =
      table.find_rows(logic::LogicalOp::And,
                      ColumnPredicate<double>{"xyzobs.px.value",
                                              [](const auto &span, size_t i) {
                                                return span(i, 0) > 200.0;
                                              }},
                      ColumnPredicate<double>{"xyzobs.px.value",
                                              [](const auto &span, size_t i) {
                                                return span(i, 2) < 5.0;
                                              }});

  std::cout << "Selected " << selected.size()
            << " rows with x > 200 AND z < 5.0\n";
  EXPECT_FALSE(selected.empty());

  // Select and validate all remaining rows satisfy both predicates
  auto filtered = table.select(selected);
  auto filtered_px = filtered.column<double>("xyzobs.px.value");
  ASSERT_TRUE(filtered_px);
  const auto &span = filtered_px.value();

  for (size_t i = 0; i < span.extent(0); ++i) {
    EXPECT_TRUE(span(i, 0) > 200.0 && span(i, 2) < 5.0);
  }
}

TEST_F(ReflectionTableTest, FindRowsThrowsOnWrongType) {
  ReflectionTable table(test_file_path.string());

  try {
    table.find_rows(
        ColumnPredicate<int>{"xyzobs.px.value", [](const auto &span, size_t i) {
                               return span(i, 0) > 0;
                             }});
    FAIL() << "Expected std::runtime_error: type mismatch";
  } catch (const std::runtime_error &e) {
    std::cout << "[FindRowsThrowsOnWrongType] Error: " << e.what() << "\n";
    EXPECT_TRUE(std::string(e.what()).find("wrong type") != std::string::npos);
  }
}
#pragma endregion

#pragma region Writing
TEST_F(ReflectionTableTest, WriteAndReloadProducesSameData) {
  // Load original table
  ReflectionTable table(test_file_path.string());

  // Take a small sample for writing
  std::vector<size_t> selected_rows = {0, 1, 2};
  ReflectionTable subset = table.select(selected_rows);

  // Temporary write path
  std::filesystem::path temp_file =
      std::filesystem::current_path() / "reflection_test_write.h5";
  std::string temp_file_str = temp_file.string();

  // Write to file
  subset.write(temp_file_str);

  // Reload written file
  ReflectionTable reloaded(temp_file_str);

  auto column_names = subset.get_column_names();
  ASSERT_EQ(column_names, reloaded.get_column_names());

  for (const auto &name : column_names) {
    // Support double and int for now
    auto original_d = subset.column<double>(name);
    auto reloaded_d = reloaded.column<double>(name);

    if (original_d && reloaded_d) {
      const auto &a = original_d.value();
      const auto &b = reloaded_d.value();

      ASSERT_EQ(a.extent(0), b.extent(0));
      ASSERT_EQ(a.extent(1), b.extent(1));

      for (size_t i = 0; i < a.extent(0); ++i) {
        for (size_t j = 0; j < a.extent(1); ++j) {
          EXPECT_DOUBLE_EQ(a(i, j), b(i, j));
        }
      }
      continue;
    }

    auto original_i = subset.column<int>(name);
    auto reloaded_i = reloaded.column<int>(name);

    if (original_i && reloaded_i) {
      const auto &a = original_i.value();
      const auto &b = reloaded_i.value();

      ASSERT_EQ(a.extent(0), b.extent(0));
      ASSERT_EQ(a.extent(1), b.extent(1));

      for (size_t i = 0; i < a.extent(0); ++i) {
        for (size_t j = 0; j < a.extent(1); ++j) {
          EXPECT_EQ(a(i, j), b(i, j));
        }
      }
      continue;
    }

    std::cerr << "Skipping unsupported type or column not found: " << name
              << "\n";
  }

  // Clean up
  std::filesystem::remove(temp_file);
}

TEST_F(ReflectionTableTest, ExperimentMetadataRoundTrip) {
  // Load original table
  ReflectionTable table(test_file_path.string());

  // Verify that metadata was loaded from the original file
  const auto &ids = table.get_experiment_ids();
  const auto &names = table.get_identifiers();

  std::cout << "\nLoaded experiment_ids:\n";
  for (const auto &id : ids) {
    std::cout << "  - " << id << "\n";
  }

  std::cout << "Loaded identifiers:\n";
  for (const auto &s : names) {
    std::cout << "  - \"" << s << "\"\n";
  }

  ASSERT_FALSE(ids.empty()) << "experiment_ids not loaded from input file.";
  ASSERT_FALSE(names.empty()) << "identifiers not loaded from input file.";
  ASSERT_EQ(ids.size(), names.size()) << "Mismatch in metadata entry count.";

  // Take a small sample with at least one experiment reference
  std::vector<size_t> selected_rows = {0, 1};
  ReflectionTable subset = table.select(selected_rows);

  // Write to temp file
  std::filesystem::path temp_file =
      std::filesystem::current_path() / "reflection_test_metadata.h5";
  std::string temp_file_str = temp_file.string();

  std::cout << "Writing subset to: " << temp_file_str << "\n";
  subset.write(temp_file_str);

  // Reload from file
  ReflectionTable reloaded(temp_file_str);

  // Compare experiment_ids
  const auto &ids_reloaded = reloaded.get_experiment_ids();
  ASSERT_EQ(ids.size(), ids_reloaded.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    std::cout << "ID[" << i << "]: original=" << ids[i]
              << " reloaded=" << ids_reloaded[i] << "\n";
    EXPECT_EQ(ids[i], ids_reloaded[i]);
  }

  // Compare identifiers
  const auto &names_reloaded = reloaded.get_identifiers();
  ASSERT_EQ(names.size(), names_reloaded.size());
  for (size_t i = 0; i < names.size(); ++i) {
    std::cout << "Identifier[" << i << "]: original=\"" << names[i]
              << "\" reloaded=\"" << names_reloaded[i] << "\"\n";
    EXPECT_EQ(names[i], names_reloaded[i]);
  }

  std::cout << "\n✅ Metadata round-trip test passed.\n";
  std::cout << "Inspect manually in HDFView if needed: " << temp_file_str
            << "\n";

  // Clean up
  std::filesystem::remove(temp_file);
}

TEST_F(ReflectionTableTest, EmptyTableWriteSucceeds) {
  ReflectionTable table;
  table.set_experiment_ids({1, 2, 3});
  table.set_identifiers({"a", "b", "c"});

  std::filesystem::path temp_file =
      std::filesystem::current_path() / "empty_table.h5";
  EXPECT_NO_THROW(table.write(temp_file.string()));

  ReflectionTable loaded(temp_file.string());
  EXPECT_TRUE(loaded.get_column_names().empty());

  std::filesystem::remove(temp_file);
}
#pragma endregion
