#include <dx2/h5/h5read_processed.h>
#include <dx2/h5/h5write.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <vector>

// Test writing a 3D vector dataset to an HDF5 file and verifying the written
// content
TEST(ExampleTests, WriteArrayTest) {
  // Define the test file path
  std::filesystem::path cwd = std::filesystem::current_path();
  std::string test_file_path = cwd.generic_string();
  test_file_path.append("/data/test_write.h5");

  // Prepare test data (3D vectors)
  std::vector<std::array<double, 3>> xyzobs_test_data = {
      {100.1, 200.2, 300.3}, {400.4, 500.5, 600.6}, {700.7, 800.8, 900.9}};

  // Call the function to write the test data to the HDF5 file
  std::string dataset_name = "xyzobs.px.value";
  std::string group_path = "processing";

  write_xyzobs_to_h5_file(test_file_path, group_path, dataset_name,
                          xyzobs_test_data);

  // Read the written data back from the HDF5 file
  std::string full_dataset_path = "/dials/processing/group_0/xyzobs.px.value";
  std::vector<double> read_xyzobs =
      read_array_from_h5_file<double>(test_file_path, full_dataset_path);

  // Validate the size of the read data
  ASSERT_EQ(read_xyzobs.size(), xyzobs_test_data.size() * 3);

  // Validate the content of the read data
  for (std::size_t i = 0; i < xyzobs_test_data.size(); ++i) {
    EXPECT_DOUBLE_EQ(read_xyzobs[i * 3 + 0], xyzobs_test_data[i][0]);
    EXPECT_DOUBLE_EQ(read_xyzobs[i * 3 + 1], xyzobs_test_data[i][1]);
    EXPECT_DOUBLE_EQ(read_xyzobs[i * 3 + 2], xyzobs_test_data[i][2]);
  }

  // Clean up test file after successful run (comment out to keep the test file)
  std::filesystem::remove(test_file_path);
}
