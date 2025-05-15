#include <array>
#include <dx2/h5/h5read_processed.hpp>
#include <dx2/h5/h5write.hpp>
#include <filesystem>
#include <gtest/gtest.h>
#include <hdf5.h>
#include <string>
#include <vector>

// Test fixture for HDF5-related tests
class HDF5Test : public ::testing::Test {
protected:
  std::filesystem::path test_file_path;

  void SetUp() override {
    test_file_path = std::filesystem::current_path() / "test_hdf5_file.h5";
  }

  void TearDown() override {
    if (std::filesystem::exists(test_file_path)) {
      std::filesystem::remove(test_file_path);
    }
  }
};

// --------------- traverse_or_create_groups TESTS ---------------
#pragma region traverse_or_create_groups tests

TEST_F(HDF5Test, TraverseOrCreateGroups) {
  // Create an HDF5 file
  hid_t file = H5Fcreate(test_file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";

  try {
    std::string group_path = "/dials/processing/group_0";
    hid_t final_group = traverse_or_create_groups(file, group_path);
    ASSERT_GE(final_group, 0) << "Failed to create or open the final group.";

    // Verify group hierarchy exists
    hid_t dials_group = H5Gopen(file, "dials", H5P_DEFAULT);
    ASSERT_GE(dials_group, 0) << "Failed to open the 'dials' group.";

    hid_t processing_group = H5Gopen(dials_group, "processing", H5P_DEFAULT);
    ASSERT_GE(processing_group, 0) << "Failed to open the 'processing' group.";

    hid_t group_0 = H5Gopen(processing_group, "group_0", H5P_DEFAULT);
    ASSERT_GE(group_0, 0) << "Failed to open the 'group_0' group.";

    // Close all groups
    H5Gclose(group_0);
    H5Gclose(processing_group);
    H5Gclose(dials_group);
    H5Gclose(final_group);
  } catch (const std::runtime_error &e) {
    FAIL() << "Runtime error occurred: " << e.what();
  }

  // Close and validate file existence
  H5Fclose(file);
  ASSERT_TRUE(std::filesystem::exists(test_file_path))
      << "HDF5 file was not created.";
}

#pragma endregion traverse_or_create_groups tests
// --------------- deduce_shape TESTS ---------------
#pragma region deduce_shape tests

TEST(DeduceShapeTests, OneDimensionalVector) {
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
  EXPECT_EQ(deduce_shape(data), (std::vector<hsize_t>{4}));
}

TEST(DeduceShapeTests, TwoDimensionalVector) {
  std::vector<std::vector<double>> data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  EXPECT_EQ(deduce_shape(data), (std::vector<hsize_t>{2, 3}));
}

TEST(DeduceShapeTests, TwoDimensionalVectorArray) {
  std::vector<std::array<double, 3>> data = {
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  EXPECT_EQ(deduce_shape(data), (std::vector<hsize_t>{3, 3}));
}

TEST(DeduceShapeTests, ThreeDimensionalVector) {
  std::vector<std::vector<std::vector<int>>> data = {{{1, 2}, {3, 4}},
                                                     {{5, 6}, {7, 8}}};
  EXPECT_EQ(deduce_shape(data), (std::vector<hsize_t>{2, 2, 2}));
}

TEST(DeduceShapeTests, EmptyContainer) {
  std::vector<std::vector<double>> data = {};
  EXPECT_EQ(deduce_shape(data), (std::vector<hsize_t>{0}));
}

TEST(DeduceShapeTests, MixedSizeInnerContainers) {
  std::vector<std::vector<double>> data = {{1.0, 2.0}, {3.0, 4.0, 5.0}};
  EXPECT_THROW(deduce_shape(data), std::exception);
}

#pragma endregion deduce_shape tests
// --------------- flatten TESTS ---------------
#pragma region flatten tests

TEST(FlattenTests, OneDimensionalVector) {
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
  EXPECT_EQ(flatten(data), data);
}

TEST(FlattenTests, TwoDimensionalVector) {
  std::vector<std::vector<double>> data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  EXPECT_EQ(flatten(data), (std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}));
}

TEST(FlattenTests, TwoDimensionalVectorArray) {
  std::vector<std::array<double, 3>> data = {
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  EXPECT_EQ(flatten(data),
            (std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}));
}

TEST(FlattenTests, ThreeDimensionalVector) {
  std::vector<std::vector<std::vector<int>>> data = {{{1, 2}, {3, 4}},
                                                     {{5, 6}, {7, 8}}};
  EXPECT_EQ(flatten(data), (std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8}));
}

TEST(FlattenTests, EmptyContainer) {
  std::vector<std::vector<double>> data = {};
  EXPECT_EQ(flatten(data), (std::vector<double>{}));
}

#pragma endregion flatten tests
// --------------- write_h5 TESTS ---------------
#pragma region write_h5 tests

TEST_F(HDF5Test, WriteOneDimensionalVector) {
  std::string dataset_path = "/group_1/dataset_1d";
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};

  write_data_to_h5_file(test_file_path.string(), dataset_path, data);
  EXPECT_EQ(
      read_array_from_h5_file<double>(test_file_path.c_str(), dataset_path),
      data);
}

TEST_F(HDF5Test, WriteTwoDimensionalVector) {
  std::string dataset_path = "/group_2/dataset_2d";
  std::vector<std::vector<double>> data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

  write_data_to_h5_file(test_file_path.string(), dataset_path, data);
  EXPECT_EQ(
      read_array_from_h5_file<double>(test_file_path.c_str(), dataset_path),
      flatten(data));
}

TEST_F(HDF5Test, WriteEmptyDataset) {
  std::string dataset_path = "/group_empty/dataset_empty";
  std::vector<double> data = {};

  write_data_to_h5_file(test_file_path.string(), dataset_path, data);
  EXPECT_TRUE(
      read_array_from_h5_file<double>(test_file_path.c_str(), dataset_path)
          .empty());
}

// Test writing to a file that already exists
TEST_F(HDF5Test, WriteToExistingFile) {
  hid_t file = H5Fcreate(test_file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";
  H5Fclose(file);

  std::string dataset_path = "/existing_group/existing_dataset";
  std::vector<double> data = {1.0, 2.0, 3.0, 4.0};

  write_data_to_h5_file(test_file_path.string(), dataset_path, data);
  EXPECT_EQ(
      read_array_from_h5_file<double>(test_file_path.c_str(), dataset_path),
      data);
}

// Test writing to a group that already exists
TEST_F(HDF5Test, WriteToExistingGroup) {
  hid_t file = H5Fcreate(test_file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";

  hid_t group =
      H5Gcreate(file, "/group_1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(group, 0) << "Failed to create group.";
  H5Gclose(group);
  H5Fclose(file);

  std::string dataset_path = "/group_1/dataset_1";
  std::vector<double> data = {10.0, 20.0, 30.0, 40.0};

  write_data_to_h5_file(test_file_path.string(), dataset_path, data);
  EXPECT_EQ(
      read_array_from_h5_file<double>(test_file_path.c_str(), dataset_path),
      data);
}

// Test writing to a dataset that already exists
TEST_F(HDF5Test, WriteToExistingDataset) {
  hid_t file = H5Fcreate(test_file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
  ASSERT_GE(file, 0) << "Failed to create HDF5 file.";

  hid_t group =
      H5Gcreate(file, "/group_2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_GE(group, 0) << "Failed to create group.";

  hsize_t dims[1] = {4};
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
  write_data_to_h5_file(test_file_path.string(), "/group_2/dataset_existing",
                        new_data);
  EXPECT_EQ(read_array_from_h5_file<double>(test_file_path.c_str(),
                                            "/group_2/dataset_existing"),
            new_data);
}

#pragma endregion write_h5 tests
