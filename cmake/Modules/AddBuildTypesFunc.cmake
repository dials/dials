# Function included by CoverageBuildConfiguration.cmake and SetDefaultBuildRelWithDebInfo.cmake

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
