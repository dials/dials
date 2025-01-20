# dx2 - Clean C++ implementation of dxtbx data model

## Quick Start

- Ensure you have an environment with CMake and compilers available (see requirements.txt)
- Dependencies:
    - HDF5 (Required)
    - gemmi (Required)
    - [Eigen](https://eigen.tuxfamily.org/) (Optional, downloaded if missing)
    - [nlohmann/json](https://github.com/nlohmann/json) (Optional, downloaded if missing)
    - [{fmt}](https://github.com/fmtlib/fmt) (Optional, downloaded if missing)
    - [GoogleTest](https://github.com/google/googletest) (Optional, downloaded if missing)

```
mkdir build
cd build
cmake ..
```
You can then build with `make`, and run the test with `./tests/test_make_models`.

## Writing And Running Tests

See https://google.github.io/googletest/primer.html. The `tests` subdirectory
has an example test.

You can run test executables directly, run `ctest`, or run `make test` (or
via your preferred build tool).

## Using dx2 from other projects

### Option 1: Already installed into the environment

In which case the normal CMake package finding should work:

```cmake
find_package(DX2)
```

### Option 2: As a subdirectory/submodule

If you want to integrate dx2 development into your current project, you can
check dx2 out into a subfolder (or submodule) of your build, and include it with:

```cmake
add_subdirectory(path/to/dx2)
```

### Option 3: Automatically use with FetchContent

If you don't need to develop dx2 but just want to use it, you can get CMake to
automatically download and include it in your project with:

```cmake
include(FetchContent)
FetchContent_Declare(
    dx2
    GIT_REPOSITORY git@github.com:dials/dx2.git
    GIT_TAG main
    FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(dx2)

...

# Example adding to a target
target_link_libraries(my-target PUBLIC dx2)
```

> [!CAUTION]
> This is linking directly to the `main` branch. This is bad practice as it
> makes it non-trivial to ensure updates, but while the project is in the
> incubation stages we don't have a public release tag to base it on.
