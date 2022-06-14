# Add compiler flags to always output coloured diagnostics
#
# This is an issue with Ninja where it captures the output to avoid
# intermixing output like the makefile generator. Unfortuantely,
# the compiler detects it's being captured and doesn't output the
# coloured diagnostics.
#
# This module adds global compile flags to turn this back on, and
# adds a configurable FORCE_COLORED_OUTPUT option.
#
# If this is called from a CMakeLists that is not the root, then it
# assumes that it's parent is being called as an add_subdirectory,
# and therefore changes to global flags are not desired,
#

# Only do this if we're being called from the root CMakeLists - if the
# parent is added as a subdirectory then it shouldn't fiddle with global
# settings
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    # Always use coloured output, unless turned off
    option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
    if (${FORCE_COLORED_OUTPUT})
        add_compile_options(
            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<C_COMPILER_ID:GNU>>:-fdiagnostics-color=always>
            $<$<OR:$<CXX_COMPILER_ID:Clang>,$<C_COMPILER_ID:Clang>>:-fcolor-diagnostics>
            $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:AppleClang>>:-fcolor-diagnostics>
            $<$<OR:$<CXX_COMPILER_ID:IntelLLVM>,$<C_COMPILER_ID:IntelLLVM>>:-fcolor-diagnostics>
        )
    endif()
endif()

