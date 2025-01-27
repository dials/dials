#include <dx2/h5/h5read_processed.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <stdlib.h>

// A test that we can read data arrays from a h5 processing file.

TEST(ExampleTests, ReadArrayTest) {

  std::string array_name = "/dials/processing/group_0/xyzobs.px.value";
  std::string flags = "/dials/processing/group_0/flags";
  // The test has been given the tests directory as the working directory
  // so we can get the path to the test file.
  std::filesystem::path cwd = std::filesystem::current_path();
  std::string fpath = cwd.generic_string();
  fpath.append("/data/cut_strong.refl");

  std::vector<double> xyzobs_px =
      read_array_from_h5_file<double>(fpath, array_name);
  // check a random value
  double expected_value = 528.86470588235295;
  EXPECT_EQ(xyzobs_px[10], expected_value);

  std::vector<std::size_t> flags_array =
      read_array_from_h5_file<std::size_t>(fpath, flags);
  // check a random value
  std::size_t expected_flag_value = 32;
  EXPECT_EQ(flags_array[5], expected_flag_value);
}
