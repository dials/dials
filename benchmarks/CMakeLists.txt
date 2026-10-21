# Google Benchmark microbenchmarks. Enable with -DDIALS_BUILD_BENCHMARKS=ON.
# Uses an installed Google Benchmark if found (e.g. conda's "benchmark"
# package), otherwise downloads it.
find_package(benchmark QUIET)
if(NOT benchmark_FOUND)
    message(STATUS "Google Benchmark not found; fetching it")
    include(FetchContent)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.9.1
    )
    FetchContent_MakeAvailable(benchmark)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(bench_connected_components bench_connected_components.cc)
target_link_libraries(bench_connected_components PUBLIC Boost::boost CCTBX::scitbx benchmark::benchmark)
