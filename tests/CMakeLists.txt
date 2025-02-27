
add_executable(test_crystal test_crystal.cxx)
target_link_libraries(test_crystal GTest::gtest_main dx2 Eigen3::Eigen nlohmann_json::nlohmann_json gemmi::gemmi_cpp)

add_executable(test_read_h5_array test_read_h5_array.cxx)
target_link_libraries(test_read_h5_array GTest::gtest_main dx2 hdf5::hdf5)

add_executable(test_write_h5_array test_write_h5_array.cxx)
target_link_libraries(test_write_h5_array GTest::gtest_main dx2 hdf5::hdf5)

include(GoogleTest)

gtest_discover_tests(test_crystal PROPERTIES LABELS dx2tests)
gtest_discover_tests(test_read_h5_array WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" PROPERTIES LABELS dx2tests)
gtest_discover_tests(test_write_h5_array WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" PROPERTIES LABELS dx2tests)
