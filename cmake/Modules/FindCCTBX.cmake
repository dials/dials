# Distributes under BSD licence

#.rst:
# FindCCTBX
# ---------
#
# Find an existing CCTBX distribution and sets up targets for using it
#
# The CCTBX distribution is found by, in order:
#
# 1. Reading the ``CCTBX_BUILD_DIR`` cache variable, if set
# 2. Reading the ``LIBTBX_BUILD`` environment variable
# 3. Using python to ``import libtbx``
#
# Components
# ^^^^^^^^^^
#
# Any ``COMPONENTS`` passed to ``find_package`` will cause this module
# to look for a libtbx-distribution module of the same name. This will
# then be exposed as the target ``CCTBX::<module_name>`` which will set
# up the include, library paths associated with that module.
#
# Other libtbx module dependencies are not currently transitively
# followed, so will need to be specified manually.
#
# Any known shared libraries for a particular tbx-module will have
# targets added themselves. For example, cctbx builds a ``libcctbx``
# shared library. If the cctbx component is requested, this will be made
# available for linking with the target name ``CCTBX::cctbx::cctbx``.
#
# The database for recognising these shared library targets can be found
# in ``module_info.json`` in the folder above this FindCCTBX module.

# If python isn't already included, pull it in here
if (NOT TARGET Python::Interpreter)
    find_package(Python COMPONENTS Interpreter REQUIRED)
endif()

# Find the location of the libtbx build directory - where libtbx_env is
function(_cctbx_determine_libtbx_build_dir)
    # Try and read it from the environment
    message(DEBUG "Looking for libtbx build dir via environment LIBTBX_BUILD")
    if (DEFINED ENV{LIBTBX_BUILD})
        get_filename_component(_ABS_BUILD "$ENV{LIBTBX_BUILD}" ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}")
        set(LIBTBX_ENV "${_ABS_BUILD}/libtbx_env")
        message(DEBUG "Checking ${LIBTBX_ENV}")
        if (EXISTS "${LIBTBX_ENV}")
            set(CCTBX_BUILD_DIR "${_ABS_BUILD}" CACHE FILEPATH "Location of CCTBX build directory")
            message(DEBUG "Found libtbx via environment LIBTBX_BUILD: ${_ABS_BUILD}")
            return()
        endif()
    endif()

    message(DEBUG "Looking for libtbx build dir via importing libtbx in python")
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import libtbx.load_env; print(abs(libtbx.env.build_path))"
                    RESULT_VARIABLE _LOAD_ENV_RESULT
                    OUTPUT_VARIABLE _LOAD_LIBTBX_BUILD_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)

    if (NOT ${_LOAD_ENV_RESULT})
        # We found it via python import
        message(DEBUG "Got libtbx build path: ${_LOAD_LIBTBX_BUILD_DIR}")
        set(CCTBX_BUILD_DIR "${_LOAD_LIBTBX_BUILD_DIR}" CACHE FILEPATH "Location of CCTBX build directory")
        return()
    endif()

    message(DEBUG "Could not find through direct python; looking for libtbx.python as last resort")
    execute_process(COMMAND "libtbx.python" -c "import libtbx.load_env; print(abs(libtbx.env.build_path))"
                    RESULT_VARIABLE _TBX_LOAD_ENV_RESULT
                    OUTPUT_VARIABLE _TBX_LOAD_LIBTBX_BUILD_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)

    if (NOT ${_TBX_LOAD_ENV_RESULT})
        # We found it via python import
        message(DEBUG "Got libtbx build path: ${_TBX_LOAD_LIBTBX_BUILD_DIR}")
        set(CCTBX_BUILD_DIR "${_TBX_LOAD_LIBTBX_BUILD_DIR}" CACHE FILEPATH "Location of CCTBX build directory")
        return()
    endif()
endfunction()

function(_read_libtbx_env RESULT_VARIABLE)
    cmake_path(SET _read_env_script NORMALIZE "${CMAKE_CURRENT_LIST_DIR}/../read_env.py")
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
        set(_read_env_is_windows "--windows")
    endif()
    # Get the system prefix from python
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys; print(sys.prefix)"
                    RESULT_VARIABLE _SYS_PREFIX_RESULT
                    OUTPUT_VARIABLE _SYS_PREFIX
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (${_SYS_PREFIX_RESULT})
        message(FATAL_ERROR "Failed to read sys.prefix out of configured python")
    endif()
    # Now, use this an read out the libtbx_env
    execute_process(
        COMMAND ${Python_EXECUTABLE}
        "${_read_env_script}"
        "${CCTBX_BUILD_DIR}/libtbx_env"
        "--build-path"
        "${CCTBX_BUILD_DIR}"
        "--sys-prefix"
        "${_SYS_PREFIX}"
        ${_read_env_is_windows}
        OUTPUT_VARIABLE _env_json
        RESULT_VARIABLE _result)
    if (_result)
        message(FATAL_ERROR "Failed to read environment file: ${CCTBX_BUILD_DIR}/libtbx_env")
    endif()
    set("${RESULT_VARIABLE}" "${_env_json}" PARENT_SCOPE)
endfunction()

# Read details for a single module out of libtbx_env and other info
function(_cctbx_read_module MODULE)
    _read_libtbx_env(_env_json)
    # We now have a json representation of libtbx_env - extract the entry for this modile
    string(JSON _module_json ERROR_VARIABLE _error GET "${_env_json}" module_dict ${MODULE})
    if (NOT _module_json)
        set("CCTBX_${MODULE}_FOUND" FALSE PARENT_SCOPE)
        return()
    else()
        set("CCTBX_${MODULE}_FOUND" TRUE PARENT_SCOPE)
    endif()
    # Now, ensure a target exists for this module
    set(_target "CCTBX::${MODULE}")
    if (NOT TARGET CCTBX::${MODULE})
        add_library(${_target} INTERFACE IMPORTED)
    endif()

    # Get the dist-wide include dir
    string(JSON _include_paths GET "${_env_json}" include_path)
    string(JSON _lib_path GET "${_env_json}" lib_path)

    # Read the metainfo database - we might have extra information we need to inject
    file(READ "${CMAKE_CURRENT_LIST_DIR}/../module_info.json" _modules_db)
    # Read a list of libraries that this module exports
    string(JSON _module_libs     ERROR_VARIABLE _error GET "${_modules_db}" "libraries" "${MODULE}")
    # Read the extra include paths for this module
    string(JSON _module_includes_array ERROR_VARIABLE _error GET "${_modules_db}" "includes" "${MODULE}")
    if (_module_includes_array)
        # Convert this array to a CMake list
        string(JSON _n_includes LENGTH "${_module_includes_array}")
        math(EXPR _n_includes "${_n_includes} - 1")  # CMake RANGE is inclusive
        foreach( _n RANGE "${_n_includes}")
            string(JSON _include GET "${_module_includes_array}" "${_n}")
            list(APPEND _module_includes "${_include}")
        endforeach()
    endif()

    # Work out what dist paths need to be consulted for this module
    string(JSON _n_dist_paths LENGTH "${_module_json}" dist_paths)
    math(EXPR _n_dist_paths "${_n_dist_paths} - 1")  # CMake RANGE is inclusive

    # We need to work out include/ directories for this module.
    #
    # Algorithm: For every listed dist path:
    #               if the module_info database contains an entry for this module:
    #                   use those
    #               else if folder has an include/ subdir:
    #                   use that
    #               else:
    #                   use the parent of the dist path above
    foreach(_n RANGE "${_n_dist_paths}")
        string(JSON _dist_path GET "${_module_json}" dist_paths ${_n})
        if(NOT _dist_path)
            continue()
        endif()
        list(APPEND _dist_paths "${_dist_path}")
        # If we don't have a specific include-path override for this,
        # then use the dist-paths as roots for the include paths
        if (_module_includes)
            # Try appending every includes dir
            foreach(_include ${_module_includes})
                # Build up a full path for this
                string(FIND "${_include}" "#build/" _build_relative)
                if ("${_build_relative}" EQUAL 0)
                    string(SUBSTRING "${_include}" 7 -1 _include)
                    # We use an include directory relative to the environment build/
                    cmake_path(APPEND CCTBX_BUILD_DIR "${_include}" OUTPUT_VARIABLE _full_include)
                else()
                    cmake_path(APPEND _dist_path "${_include}" OUTPUT_VARIABLE _full_include)
                endif()
                cmake_path(ABSOLUTE_PATH _full_include NORMALIZE)
                # We might have multiple dist paths. Only use include dirs that exist.
                if(EXISTS "${_full_include}")
                    list(APPEND _include_paths "${_full_include}")
                endif()
            endforeach()
        else()
            # We didn't have any specific override. If include/ exists
            # then use that, otherwise use the module parent directory.
            if (EXISTS "${_dist_path}/include")
                list(APPEND _include_paths "${_dist_path}/include")
            else()
                cmake_path(GET _dist_path PARENT_PATH _result)
                list(APPEND _include_paths "${_result}")
            endif()
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _include_paths)
    list(REMOVE_DUPLICATES _dist_paths)
    target_include_directories(${_target} INTERFACE ${_include_paths})
    message(DEBUG "${MODULE} include directories: ${_include_paths}")
    set(CCTBX_${MODULE}_DIST "${_dist_paths}" PARENT_SCOPE)

    # Find if this module has a "base" library, and any sub-libraries
    # string(JSON _modules_libs GET "${_modules_libs}" "libraries")
    if (_module_libs)
        # We have actual libraries to import as imported libraries
        # iterate over every key: value in the object
        string(JSON _n_libs LENGTH "${_module_libs}")
        math(EXPR _n_libs "${_n_libs} - 1")
        foreach(_n RANGE ${_n_libs})
            string(JSON _name MEMBER "${_module_libs}" ${_n})
            string(JSON _libnames GET "${_module_libs}" "${_name}")
            set(_lib_searchname "_lib_${MODULE}_${_name}")


            # Find this library - or the first of a list of fallback options
            foreach(_libname ${_libnames})
                message(DEBUG "Processing ${_libname}")
                set(lib_specific_name "_liblocation_CCTBX_${MODULE}_${_libname}")
                find_library(${lib_specific_name} ${_libname} HINTS ${_lib_path})
                if(${lib_specific_name})
                    set(${_lib_searchname} "${${lib_specific_name}}")
                    message(DEBUG "Found ${lib_specific_name}=${${lib_specific_name}}")
                    break()
                endif()
                message(DEBUG "Didn't find lib${_libname} for ${MODULE}::${_name}")
            endforeach()

            if (NOT ${_lib_searchname})
                # If this library isn't present, it might not be important - so warn only
                message(WARNING "Libtbx module ${MODULE} has library named lib${_libname} but cannot find it - the module may be misconfigured")
            else()
                if(_name STREQUAL base)
                    message(DEBUG "Module library ${_target} has root library at ${${_lib_searchname}}")
                    target_link_libraries(${_target} INTERFACE "${${_lib_searchname}}")
                else()
                    message(DEBUG "Extra library ${_target}::${_name} = ${${_lib_searchname}}")
                    add_library(${_target}::${_name} INTERFACE IMPORTED)
                    target_link_libraries(${_target}::${_name} INTERFACE ${_target} "${${_lib_searchname}}")
                endif()
            endif()
        endforeach()
    endif()
endfunction()


if (NOT CCTBX_BUILD_DIR)
    _cctbx_determine_libtbx_build_dir()
endif()
# Make this absolute, in case it was specified as relative
cmake_path(ABSOLUTE_PATH CCTBX_BUILD_DIR BASE_DIRECTORY "${CMAKE_BINARY_DIR}" NORMALIZE)

if (CCTBX_BUILD_DIR)
    message(DEBUG "Using build dir ${CCTBX_BUILD_DIR}")

    foreach(_comp IN LISTS CCTBX_FIND_COMPONENTS)
        # Only try and find this component if we didn't already find it
        if(NOT TARGET CCTBX::${_comp})
            _cctbx_read_module(${_comp})
        else()
            # We already had a target for this module, so mark found
            set("CCTBX_${_comp}_FOUND" TRUE)
        endif()
    endforeach()
    unset(_comp)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CCTBX
    REQUIRED_VARS CCTBX_BUILD_DIR
    HANDLE_COMPONENTS
)
