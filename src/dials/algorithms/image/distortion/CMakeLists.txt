Python_add_library(
    dials_algorithms_image_distortion_ext
    MODULE
    boost_python/distortion_ext.cc
)
target_link_libraries( dials_algorithms_image_distortion_ext PUBLIC CCTBX::cctbx Boost::python )