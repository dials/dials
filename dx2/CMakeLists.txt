add_library(
    dx2
    empty.cxx
)
target_link_libraries(
    dx2
    PRIVATE
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    hdf5::hdf5
    gemmi::gemmi_cpp
)

target_include_directories(
    dx2 PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(dx2 PUBLIC $<TARGET_NAME_IF_EXISTS:hdf5::hdf5>)
