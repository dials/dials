
add_executable(test_example test_example.cxx)
target_link_libraries(test_example GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_example)