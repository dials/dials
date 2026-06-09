Python_add_library(
    dials_algorithms_image_connected_components_ext
    MODULE
    boost_python/connected_components_ext.cc
)
target_link_libraries( dials_algorithms_image_connected_components_ext PUBLIC CCTBX::cctbx Boost::python )