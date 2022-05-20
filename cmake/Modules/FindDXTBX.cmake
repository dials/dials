# Distributes under BSD licence

#.rst:
# FindDXTBX
# ---------
#
# Finds the DXTBX package, via CCTBX if necessary.
#
# Creates DXTBX::dxtbx target.

# Don't reconfigure every time
if (TARGET DXTBX::DXTBX)
    return()
endif()

# If python isn't already included, pull it in here
if (NOT TARGET Python::Interpreter)
    find_package(Python COMPONENTS Interpreter REQUIRED)
endif()

message(DEBUG "Looking for dxtbx build dir via importing in python")
execute_process(COMMAND ${Python_EXECUTABLE} -c "import dxtbx, pathlib; print(pathlib.Path(dxtbx.__file__).parent.resolve())"
                RESULT_VARIABLE _IMPORT_RESULT
                OUTPUT_VARIABLE _IMPORT_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET)

if (NOT ${_IMPORT_RESULT})
    # We found it via python import
    set(DXTBX_DIR "${_IMPORT_DIR}")
else()
    message(DEBUG "Could not import dxtbx directly; Searching via CCTBX")
    # Try to find it via CCTBX module
    find_package(CCTBX QUIET REQUIRED COMPONENTS dxtbx)
    if (CCTBX_FOUND)
        if(NOT CCTBX_dxtbx_DIST)
            message(FATAL_ERROR "Could not determine dxtbx location; CCTBX returned empty dist path")
        endif()
        set(DXTBX_DIR "${CCTBX_dxtbx_DIST}/src/dxtbx")
    endif()
endif()

if (DXTBX_DIR)
    add_library(DXTBX::DXTBX INTERFACE IMPORTED)
    cmake_path(GET DXTBX_DIR PARENT_PATH _dxtbx_include_root)
    target_include_directories(DXTBX::DXTBX INTERFACE "${_dxtbx_include_root}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DXTBX
    REQUIRED_VARS DXTBX_DIR
    HANDLE_COMPONENTS
)
