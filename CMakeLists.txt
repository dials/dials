cmake_minimum_required(VERSION 3.28...3.30)

project(dx2 LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only include global build-configuration features if we are the root CMakeLists.txt
if(CMAKE_CURRENT_SOURCE_DIR EQUAL CMAKE_SOURCE_DIR)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
    include(SetDefaultBuildRelWithDebInfo)
    include(AlwaysColourCompilation)
endif()

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# #######################################################################
# External Dependencies
find_package(HDF5 REQUIRED)
find_package(gemmi CONFIG REQUIRED)

# #######################################################################
# Automatic Dependencies
set(FETCHCONTENT_QUIET OFF)
include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
    EXCLUDE_FROM_ALL
    FIND_PACKAGE_ARGS
)
FetchContent_Declare(
    GTest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG b514bdc898e2951020cbdca1304b75f5950d1f59 # v1.15.2
    EXCLUDE_FROM_ALL
    FIND_PACKAGE_ARGS
)
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 11.0.2
    EXCLUDE_FROM_ALL
    FIND_PACKAGE_ARGS
)
FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    EXCLUDE_FROM_ALL
    FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(Eigen3 fmt nlohmann_json)

# GTest could have been made available under a different name
if(NOT TARGET GTest::gtest_main)
    FetchContent_MakeAvailable(GTest)
endif()

# #######################################################################
# Targets
enable_testing()

add_subdirectory(dx2)
add_subdirectory(tests)

# #######################################################################
# Installation
install(
    TARGETS dx2
    EXPORT DX2Targets
    LIBRARY DESTINATION lib

    # FILE_SET HEADERS DESTINATION include/dx2
)
install(
    FILES
    include/dx2/dx2.h
    DESTINATION include/dx2
)

install(
    EXPORT DX2Targets
    FILE DX2Targets.cmake
    DESTINATION lib/cmake/DX2
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/DX2Config.cmake"
    INSTALL_DESTINATION "lib/cmake/DX2"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DX2ConfigVersion.cmake"
    VERSION "${DX2_VERSION_MAJOR}.${DX2_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
)
install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/DX2Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DX2ConfigVersion.cmake
    DESTINATION lib/cmake/DX2
)
export(EXPORT DX2Targets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/DX2Targets.cmake"
)
