#include <dx2/h5/h5write.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <hdf5.h>
#include <string>

// Test the traverse_or_create_groups function
TEST(HDF5Tests, TraverseOrCreateGroupsTest) {
    // Define the test file path
    std::filesystem::path cwd = std::filesystem::current_path();
    std::string test_file_path = cwd.generic_string();
    test_file_path.append("/data/test_traverse_or_create_groups.h5");

    // Create or open an HDF5 file
    hid_t file = H5Fcreate(test_file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
