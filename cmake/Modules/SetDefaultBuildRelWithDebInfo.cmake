# If unset, set the default build type to RelWithDebInfo
#
# This also adds a proper set of property strings to the
# CMAKE_BUILD_TYPE configuration variable.
#
# Nothing will be done if the including script is not the root script.
#
# From https://blog.kitware.com/cmake-and-the-default-build-type/

# Only do this if we're being called from the root CMakeLists - if the
# parent is added as a subdirectory then it shouldn't fiddle with global
# settings
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  # Set the default build type to RelWithDebInfo
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()
endif()
