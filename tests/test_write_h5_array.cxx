#include <array>
#include <dx2/h5/h5read_processed.h>
#include <dx2/h5/h5write.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <hdf5.h>
#include <string>
#include <vector>

#pragma region traverse_or_create_groups tests
// Test the traverse_or_create_groups function
TEST(HDF5Tests, TraverseOrCreateGroupsTest) {
  // Define the test file path
  std::filesystem::path cwd = std::filesystem::current_path();
  std::string test_file_path = cwd.generic_string();
  test_file_path.append("/data/test_traverse_or_create_groups.h5");

  // Create or open an HDF5 file
  hid_t file = H5Fcreate(test_file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";

  try {
    // Define a hierarchical path
    std::string group_path = "/dials/processing/group_0";

    // Call the function to traverse or create groups
    hid_t final_group = traverse_or_create_groups(file, group_path);
    ASSERT_GE(final_group, 0) << "Failed to create or open the final group.";

    // Verify that each group in the hierarchy exists
    hid_t dials_group = H5Gopen(file, "dials", H5P_DEFAULT);
    ASSERT_GE(dials_group, 0) << "Failed to open the 'dials' group.";

    hid_t processing_group = H5Gopen(dials_group, "processing", H5P_DEFAULT);
    ASSERT_GE(processing_group, 0) << "Failed to open the 'processing' group.";

    hid_t group_0 = H5Gopen(processing_group, "group_0", H5P_DEFAULT);
    ASSERT_GE(group_0, 0) << "Failed to open the 'group_0' group.";

    // Close all opened groups
    H5Gclose(group_0);
    H5Gclose(processing_group);
    H5Gclose(dials_group);
    H5Gclose(final_group);
  } catch (const std::runtime_error &e) {
    FAIL() << "Runtime error occurred: " << e.what();
  }

  // Close the file
  H5Fclose(file);

  // Validate that the HDF5 file was created
  ASSERT_TRUE(std::filesystem::exists(test_file_path))
      << "HDF5 file was not created.";

  // Clean up test file after successful run (comment out to keep the test file)
  std::filesystem::remove(test_file_path);
}
#pragma endregion traverse_or_create_groups tests

#pragma region deduce_shape tests
// Test deducing shape for 1D container (std::vector)
TEST(DeduceShapeTests, OneDimensionalVector) {
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
  std::vector<hsize_t> expected_shape = {4}; // Single dimension of size 4
  auto shape = deduce_shape(data);
  EXPECT_EQ(shape, expected_shape);
}

// Test deducing shape for 2D container (std::vector<std::vector>)
TEST(DeduceShapeTests, TwoDimensionalVector) {
  std::vector<std::vector<double>> data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  std::vector<hsize_t> expected_shape = {2, 3}; // 2 rows, 3 columns
  auto shape = deduce_shape(data);
  EXPECT_EQ(shape, expected_shape);
}

// Test deducing shape for 2D container (std::vector<std::array>)
TEST(DeduceShapeTests, TwoDimensionalVectorArray) {
  std::vector<std::array<double, 3>> data = {
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  std::vector<hsize_t> expected_shape = {3, 3}; // 3 rows, 3 columns
  auto shape = deduce_shape(data);
  EXPECT_EQ(shape, expected_shape);
}

// Test deducing shape for 3D container (std::vector<std::vector<std::vector>>)
TEST(DeduceShapeTests, ThreeDimensionalVector) {
  std::vector<std::vector<std::vector<int>>> data = {{{1, 2}, {3, 4}},
                                                     {{5, 6}, {7, 8}}};
  std::vector<hsize_t> expected_shape = {2, 2, 2}; // 2x2x2 structure
  auto shape = deduce_shape(data);
  EXPECT_EQ(shape, expected_shape);
}

// Test deducing shape for empty container
TEST(DeduceShapeTests, EmptyContainer) {
  std::vector<std::vector<double>> data = {};
  std::vector<hsize_t> expected_shape = {0}; // Outer container size is 0
  auto shape = deduce_shape(data);
  EXPECT_EQ(shape, expected_shape);
}

// Test deducing shape for mixed-sized inner containers (should throw or handle
// gracefully)
TEST(DeduceShapeTests, MixedSizeInnerContainers) {
  std::vector<std::vector<double>> data = {
      {1.0, 2.0}, {3.0, 4.0, 5.0} // Different size
  };

  // deduce_shape assumes uniformity. If it encounters mixed sizes, it might
  // throw or return the first subcontainer's size.
  try {
    auto shape = deduce_shape(data);
    FAIL() << "Expected deduce_shape to throw due to inconsistent inner "
              "container sizes.";
  } catch (const std::exception &e) {
    SUCCEED(); // deduce_shape threw an exception as expected
  } catch (...) {
    FAIL() << "Expected deduce_shape to throw a standard exception.";
  }
}
#pragma endregion deduce_shape tests

#pragma region flatten tests
// Test flattening a 1D vector
TEST(FlattenTests, OneDimensionalVector) {
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
  auto flat_data = flatten(data); // Should return the same as the input

  std::vector<double> expected_flat_data = {1.0, 2.0, 3.0, 4.0};
  EXPECT_EQ(flat_data, expected_flat_data);
}

// Test flattening a 2D vector
TEST(FlattenTests, TwoDimensionalVector) {
  std::vector<std::vector<double>> data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  auto flat_data = flatten(data);

  std::vector<double> expected_flat_data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  EXPECT_EQ(flat_data, expected_flat_data);
}

// Test flattening a 2D vector of arrays
TEST(FlattenTests, TwoDimensionalVectorArray) {
  std::vector<std::array<double, 3>> data = {
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  auto flat_data = flatten(data);

  std::vector<double> expected_flat_data = {1.0, 2.0, 3.0, 4.0, 5.0,
                                            6.0, 7.0, 8.0, 9.0};
  EXPECT_EQ(flat_data, expected_flat_data);
}

// Test flattening a 3D vector
TEST(FlattenTests, ThreeDimensionalVector) {
  std::vector<std::vector<std::vector<int>>> data = {{{1, 2}, {3, 4}},
                                                     {{5, 6}, {7, 8}}};
  auto flat_data = flatten(data);

  std::vector<int> expected_flat_data = {1, 2, 3, 4, 5, 6, 7, 8};
  EXPECT_EQ(flat_data, expected_flat_data);
}

// Test flattening an empty container
TEST(FlattenTests, EmptyContainer) {
  std::vector<std::vector<double>> data = {};
  auto flat_data = flatten(data);

  std::vector<double> expected_flat_data = {};
  EXPECT_EQ(flat_data, expected_flat_data);
}

// Test flattening a container with mixed sizes (shouldn't be an issue here)
TEST(FlattenTests, MixedSizeInnerContainers) {
  std::vector<std::vector<double>> data = {{1.0, 2.0}, {3.0, 4.0, 5.0}};

  // Flatten should concatenate all elements into a single 1D vector
  auto flat_data = flatten(data);
  std::vector<double> expected_flat_data = {1.0, 2.0, 3.0, 4.0, 5.0};
  EXPECT_EQ(flat_data, expected_flat_data);
}
#pragma endregion flatten tests

#pragma region write_h5 tests
// Test writing a 1D vector to an HDF5 file
TEST(WriteDataTests, WriteOneDimensionalVector) {
  std::string filename = "test_1d_vector.h5";
  std::string dataset_path = "/group_1/dataset_1d";

  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
  write_data_to_h5_file(filename, dataset_path, data);

  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);

  EXPECT_EQ(data, read_data);

  // Clean up test file
  std::filesystem::remove(filename);
}

// Test writing a 2D vector to an HDF5 file
TEST(WriteDataTests, WriteTwoDimensionalVector) {
  std::string filename = "test_2d_vector.h5";
  std::string dataset_path = "/group_2/dataset_2d";

  std::vector<std::vector<double>> data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  write_data_to_h5_file(filename, dataset_path, data);

  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);

  // Flatten the original data for comparison
  std::vector<double> expected_data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  EXPECT_EQ(expected_data, read_data);

  // Clean up test file
  std::filesystem::remove(filename);
}

// Test writing a 2D vector of arrays to an HDF5 file
TEST(WriteDataTests, WriteTwoDimensionalVectorArray) {
  std::string filename = "test_2d_vector_array.h5";
  std::string dataset_path = "/group_3/dataset_2d_array";

  std::vector<std::array<double, 3>> data = {
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  write_data_to_h5_file(filename, dataset_path, data);

  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);

  // Flatten the original data for comparison
  std::vector<double> expected_data = {1.0, 2.0, 3.0, 4.0, 5.0,
                                       6.0, 7.0, 8.0, 9.0};
  EXPECT_EQ(expected_data, read_data);

  // Clean up test file
  std::filesystem::remove(filename);
}

// Test writing an empty dataset to an HDF5 file
TEST(WriteDataTests, WriteEmptyDataset) {
  std::string filename = "test_empty_dataset.h5";
  std::string dataset_path = "/group_empty/dataset_empty";

  std::vector<double> data = {}; // Empty dataset
  write_data_to_h5_file(filename, dataset_path, data);

  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);

  EXPECT_TRUE(read_data.empty());

  // Clean up test file
  std::filesystem::remove(filename);
}

// Test writing to a file that already exists
TEST(WriteDataTests, WriteToExistingFile) {
  std::string filename = "test_existing_file.h5";
  std::string dataset_path = "/existing_group/existing_dataset";

  // Create the file first
  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";
  H5Fclose(file);

  // Now attempt to write data into it
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
  write_data_to_h5_file(filename, dataset_path, data);

  // Read back the data to verify it was written correctly
  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);
  EXPECT_EQ(data, read_data);

  // Clean up test file
  std::filesystem::remove(filename);
}

// Test writing to a group that already exists
TEST(WriteDataTests, WriteToExistingGroup) {
  std::string filename = "test_existing_group.h5";
  std::string dataset_path = "/group_1/dataset_1";

  // Create file and group manually
  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";

  hid_t group =
      H5Gcreate(file, "/group_1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(group, 0) << "Failed to create group.";

  H5Gclose(group);
  H5Fclose(file);

  // Now attempt to write data into the existing group
  std::vector<double> data = {10.0, 20.0, 30.0, 40.0};
  write_data_to_h5_file(filename, dataset_path, data);

  // Read back the data to verify it was written correctly
  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);
  EXPECT_EQ(data, read_data);

  // Clean up test file
  std::filesystem::remove(filename);
}

// Test writing to a dataset that already exists
TEST(WriteDataTests, WriteToExistingDataset) {
  std::string filename = "test_existing_dataset.h5";
  std::string dataset_path = "/group_2/dataset_existing";

  // Create file and dataset manually
  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";

  hid_t group =
      H5Gcreate(file, "/group_2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(group, 0) << "Failed to create group.";

  hsize_t dims[1] = {4}; // 1D dataset with 4 elements
  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset = H5Dcreate(group, "dataset_existing", H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(dataset, 0) << "Failed to create dataset.";

  std::vector<double> initial_data = {1.1, 2.2, 3.3, 4.4};
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           initial_data.data());

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Gclose(group);
  H5Fclose(file);

  // Now attempt to overwrite the dataset
  std::vector<double> new_data = {5.5, 6.6, 7.7, 8.8};
  write_data_to_h5_file(filename, dataset_path, new_data);

  // Read back the data to verify it was overwritten correctly
  auto read_data = read_array_from_h5_file<double>(filename, dataset_path);
  EXPECT_EQ(new_data, read_data);

  // Clean up test file
  std::filesystem::remove(filename);
}
#pragma endregion write_h5 tests