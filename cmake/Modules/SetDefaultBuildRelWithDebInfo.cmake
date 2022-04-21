# If unset, set the default build type to RelWithDebInfo
#
# This also adds a proper set of property strings to the
# CMAKE_BUILD_TYPE configuration variable.
#
# Nothing will be done if the including script is not the root script.
#
# From https://blog.kitware.com/cmake-and-the-default-build-type/

# Add a set of strings to the build_type metadata description, if not present
function(_add_build_types_to_cache_strings)
    get_property(_config_list CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS)
    foreach(build_type IN ITEMS "${ARGN}")
        if (NOT "${build_type}" IN_LIST _config_list)
            list(APPEND _config_list "${build_type}")
        endif()
    endforeach()
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${_config_list}")
endfunction()

# Only do this if we're being called from the root CMakeLists - if the
# parent is added as a subdirectory then it shouldn't fiddle with global
# settings
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    # Set the default build type to RelWithDebInfo
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
        _add_build_types_to_cache_strings(Debug Release MinSizeRel RelWithDebInfo)
    endif()
endif()
